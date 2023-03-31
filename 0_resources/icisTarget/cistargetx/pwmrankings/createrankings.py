#!/usr/bin/env python
import ConfigParser
import argparse
import math
import multiprocessing
import operator
import os
import re
import sys
import uuid
import zipfile

from io import BytesIO

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

import orderstatistics
from cistargetx.common.progressmonitor import ProgressMonitor


# Constants.
PROGRESS_FILE_PREFIX = "processed-motifs-"
PROGRESS_FILE_EXTENSION = ".lst"


def cbust_relative_rankings_iterator(cfg, regions, motif_filename, fasta_filename, species2region2max_score=None,
                                     cur_species=None):
    """ Iterator over (region, relative_ranking) tuples. """

    # Running Cluster Buster.
    exec_command = [cfg.get('cbust', 'command'), '-f3', '-c0',
                    os.path.join(cfg.get('cbust', 'motif_folder'), motif_filename),
                    os.path.join(cfg.get('regions', 'fasta_folder'), fasta_filename)]

    try:
        stdout, stderr = subprocess.Popen(exec_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    except OSError, msg:
        # CAVE: do not check stderr, because this stream is used by Cluster Buster for outputting progress information.
        print >> sys.stderr, "cbust: error during execution of '" + ' '.join(exec_command) + "': " + str(msg)
        sys.exit(1)

    # Parsing output of Cluster Buster:
    #
    #    cluster score    cluster start    cluster stop    region id                       motif score
    #    8.06e-06         78753            78768           intron1#NM_006042(HS3ST3A1)     0.451
    #    10.4             80292            80429           intron1#NM_006042(HS3ST3A1)     14.3
    #
    #   >>> match.groups()
    #   ('8.06e-06', 'e-06', '78755', '78768', 'intron1#NM_006042(HS3ST3A1)')
    #
    #   >>> match.groups()
    #   ('10.4', None, '80292', '80429', 'intron1#NM_006042(HS3ST3A1)')

    pattern = re.compile(r"([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+(\d+)\s+(\d+)\s+(\S+)\s+.*")
    region2max_score = dict(zip(regions, [0.0] * len(regions)))

    with BytesIO(stdout) as clusterbuster_fh:
        for line in clusterbuster_fh:
            match = pattern.match(line)

            # Skip lines without scores.
            if not match:
                continue

            # Using score of homotypic cluster.
            score = float(match.group(1))
            region = match.group(5)
            region2max_score[region] = max(region2max_score[region], score)

    if cur_species and isinstance(species2region2max_score, dict):
        species2region2max_score[cur_species] = region2max_score

    # Creating relative order of regions.
    #
    # CAVE: Some regions have the same score, these regions are ranked consecutively but the exact order is not fixed.
    # CAVE: Comparison with 0.0 is avoided. Instead greater than 0.0 is used,
    #       which works because a Cluster Buster score is never negative.
    scores = sorted(region2max_score.items(), reverse=True, key=operator.itemgetter(1))
    regions_count = float(len(scores))
    zero_score_count = float(sum(1 for t in scores if not t[1] > 0.0))
    bottom_ranking = (regions_count - (zero_score_count / 2.0)) / regions_count

    for idx, t in enumerate(scores):
        region, score = t
        rank = float(idx + 1) / regions_count if score > 0.0 else bottom_ranking

        yield region, rank


def log(cfg, msg, log_file=None):
    verbose = cfg.getboolean('log', 'verbose')
    multiprocessing_enabled = cfg.getint('multiprocessing', 'number_of_cpus') > 1
    msg = "process-ID = " + str(os.getpid()) + "; " + msg if multiprocessing_enabled else msg

    if verbose:
        print msg
        sys.stdout.flush()

    if log_file:
        print >> log_file, msg
        log_file.flush()


def save_temporary_files(cfg, pwm, output_dir, archives_filenames):
    if not cfg.getboolean('cbust', 'intermediate_files'):
        for filename in archives_filenames:
            os.remove(filename)
    else:
        zip_filename = os.path.join(output_dir, pwm + ".zip")
        zip_output_fh = zipfile.ZipFile(zip_filename, mode='w', compression=zipfile.ZIP_DEFLATED, allowZip64=True)

        for filename in archives_filenames:
            zip_output_fh.write(filename, arcname=os.path.split(filename)[1])
            os.remove(filename)

        zip_output_fh.close()


def load_gene2regions(filename):
    """ Load lookup table that maps a gene ID to its associated region IDs. """
    gene2regions = dict()

    with open(filename, 'r') as input_fh:
        for line in input_fh:
            columns = line.rstrip().split('\t')
            gene = columns[0]

            for region in columns[1:]:
                gene2regions.setdefault(gene, set()).add(region)

    return gene2regions


def write_regions_score_table(filename, species2region2max_score):
    species = sorted(species2region2max_score.keys())

    if not species:
        return

    regions = sorted(species2region2max_score[species[0]].keys())

    with open(filename, 'w') as output_fh:
        output_fh.write('.\t' + '\t'.join(species) + '\n')

        for cur_region in regions:
            output_fh.write(cur_region + '\t' + '\t'.join([str(species2region2max_score[cur_species][cur_region])
                                                           for cur_species in species]) + '\n')


def compute_rankings(cfg, pwms, species, regions, output_dir, log_file, progress_file):
    gene2regions = (load_gene2regions(cfg.get('regions', 'lut_file'))
                    if cfg.has_option('regions', 'lut_file') and cfg.get('regions', 'lut_file')
                    else None)
    os.nice(cfg.getint('multiprocessing', 'nice_increment'))
    pm = ProgressMonitor(len(pwms) * (len(species) + 1), len(species) + 1)

    for pwm, motif_filename in pwms:
        ranking_filenames = []
        archive_filenames = []
        species2region2max_score = dict() if cfg.getboolean('cbust', 'intermediate_files') else None

        for cur_species, fasta_filename in species:
            pm.start_iteration()

            log(cfg, "running cbust to score " + cur_species + " for " + pwm + " " + pm.progress_str, log_file)

            # Running Cluster Buster to score FASTA file for homotypic cluster of motif
            # and saving relative order to file.
            filename = os.path.join(output_dir, pwm + "." + cur_species + ".rr")
            ranking_filenames.append(filename)
            archive_filenames.append(filename)

            with open(filename, 'w') as output_fh:
                for region, rank in sorted(cbust_relative_rankings_iterator(cfg,
                                                                            regions,
                                                                            motif_filename,
                                                                            fasta_filename,
                                                                            species2region2max_score,
                                                                            cur_species),
                                           reverse=False,
                                           key=operator.itemgetter(1)):
                    # CAVE: Rank must be saved using g format specifier.
                    print >> output_fh, "{0:s}\t{1:.16g}".format(region, rank)

            pm.end_iteration()

        if cfg.getboolean('cbust', 'intermediate_files'):
            table_filename = os.path.join(output_dir, pwm + ".scores")
            write_regions_score_table(table_filename, species2region2max_score)
            archive_filenames.append(table_filename)

        pm.start_iteration()

        log(cfg, "running orderstatistics for " + pwm + " " + pm.progress_str, log_file)

        # Running order statistics on rankings and saving relative order to file.
        filename = os.path.join(output_dir, pwm + ".r")

        with open(filename, 'w') as output_fh:
            for feature_id, rank in sorted(orderstatistics.combined_rankings_iterator(cfg,
                                                                                      ranking_filenames,
                                                                                      gene2regions),
                                           reverse=False,
                                           key=operator.itemgetter(0)):
                print >> output_fh, "{0:s}\t{1:d}".format(feature_id, rank)

        # Deleting intermediate files if necessary.
        save_temporary_files(cfg, pwm, output_dir, archive_filenames)
        print >> progress_file, pwm
        progress_file.flush()
        pm.end_iteration()

    if log_file:
        log_file.close()

    progress_file.close()


def get_processed_motifs(output_dir):
    result = set()

    for name in os.listdir(output_dir):
        if name.startswith(PROGRESS_FILE_PREFIX) and name.endswith(PROGRESS_FILE_EXTENSION):
            with open(os.path.join(output_dir, name), 'r') as input_fh:
                for line in input_fh.readlines():
                    result.add(line.rstrip())

    return result


def main():
    parser = argparse.ArgumentParser(prog='ctx-pwm2r',
                                     description='Score Cluster-Buster PWMs across multiple species to create rankings files.')

    parser.add_argument('-i', '--inifile',
                        required=True,
                        action='store',
                        type=str,
                        dest='ini_file',
                        help='config file')
    parser.add_argument('-o', '--outputdir',
                        required=True,
                        action='store',
                        type=str,
                        dest='output_dir',
                        help='output directory which will contain the rankings files')
    parser.add_argument('-m', '--motiflistfile',
                        required=False,
                        action='store',
                        type=str,
                        dest='motifs_list_file',
                        help='Restrict list of PWMs to score from PMWs available in the motif directory to the list of motif IDs provided in this file')

    args = parser.parse_args()

    ini_file = args.ini_file

    if not os.path.isfile(ini_file):
        print >> sys.stderr, 'ERROR: Config file "{0:s}" does not exist.'.format(ini_file)
        sys.exit(1)

    cfg = ConfigParser.RawConfigParser()
    cfg.read(ini_file)

    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    processed_pwms = get_processed_motifs(output_dir)
    processed_pwms_count = len(processed_pwms)

    cbust_motif_extension = '.' + cfg.get('cbust', 'motif_extension')
    cbust_motif_extension_length = len(cbust_motif_extension)

    def get_pwm_name(name):
        return name[0:-cbust_motif_extension_length]

    requested_pwms = set()

    if args.motifs_list_file:
        if not os.path.isfile(args.motifs_list_file):
            print >> sys.stderr, 'ERROR: Motif list file "{0:s}" does not exist.'.format(args.motifs_list_file)
            sys.exit(1)

        with open(args.motifs_list_file, 'r') as fh:
            for line in fh:
                motif_id = line.strip()

                if motif_id == '' or motif_id.startswith('#'):
                    continue

                if motif_id.endswith(cbust_motif_extension):
                    requested_pwms.add(get_pwm_name(motif_id))
                else:
                    requested_pwms.add(motif_id)

    if requested_pwms:
        # Only process requested PWMs from motif list file (if they are not processed already).
        pwms = [(get_pwm_name(name), name)
                for name in os.listdir(cfg.get('cbust', 'motif_folder'))
                if (name.endswith(cbust_motif_extension) and
                    get_pwm_name(name) not in processed_pwms and
                    get_pwm_name(name) in requested_pwms)]
    else:
        # Process all PWMs (if they are not processed already) if no motif list file is used.
        pwms = [(get_pwm_name(name), name)
                for name in os.listdir(cfg.get('cbust', 'motif_folder'))
                if (name.endswith(cbust_motif_extension) and
                    get_pwm_name(name) not in processed_pwms)]

    species = [(name[0:-(1 + len(cfg.get('regions', 'fasta_extension')))], name)
               for name in os.listdir(cfg.get('regions', 'fasta_folder'))
               if name.endswith("." + cfg.get('regions', 'fasta_extension'))]

    assert len(species) > 0, "The name of the FASTA files should be encoded using the following format: <species>-<delineation>.<extension>"

    with open(cfg.get('regions', 'region_ids_file'), 'r') as input_fh:
        regions = map(str.rstrip, input_fh.readlines())

    if processed_pwms_count > 0:
        remaining_pwms_count = len(pwms)
        log(cfg,
            "{0:d}# pwms are already processed. Calculating rankings for remaining {1:d}# pwms".format(
                processed_pwms_count, remaining_pwms_count)
        )

    # One of the imported python modules sets the CPU affinity to one core.
    # If we do not recent the CPU affinity all processes (python and cbust) wil run
    # on the same core, defeating the whole purpose of the multiprocessing module.
    #
    # Most likely one of the imported modules in linked with OpenBLAS:
    #     http://stackoverflow.com/questions/23537716/importing-scipy-breaks-multiprocessing-support-in-python
    #
    # If we ever upgrade to Python 3, the call to taskset can be replaced with:
    #     os.sched_setaffinity(pid, mask)
    DEVNULL = open(os.devnull, 'w')
    subprocess.call(['taskset','-p', '0xFFFFFFFF', str(os.getpid())], stdout=DEVNULL)
    DEVNULL.close()

    process_nr = cfg.getint('multiprocessing', 'number_of_cpus')
    log_filename_extension = cfg.get('log', 'output_file_extension')
    log_filename_prefix = cfg.get('log', 'output_file_prefix')
    processes = []
    block_size = int(math.ceil(float(len(pwms)) / float(process_nr)))

    for idx in range(process_nr):
        log_output_fh = open(os.path.join(output_dir,
                                          log_filename_prefix + "-" + str(idx) + "." + log_filename_extension),
                             'w')
        progress_output_fh = open(os.path.join(output_dir,
                                               PROGRESS_FILE_PREFIX + str(uuid.uuid1()) + PROGRESS_FILE_EXTENSION),
                                  'w')

        start_idx = idx * block_size
        end_idx = min((idx + 1) * block_size, len(pwms))
        pwm_partition = pwms[start_idx:end_idx]

        p = multiprocessing.Process(target=compute_rankings,
                                    args=(cfg, pwm_partition, species, regions, output_dir, log_output_fh,
                                          progress_output_fh))
        p.start()
        processes.append(p)

    for idx in range(process_nr):
        processes[idx].join()


if __name__ == "__main__":
    main()
