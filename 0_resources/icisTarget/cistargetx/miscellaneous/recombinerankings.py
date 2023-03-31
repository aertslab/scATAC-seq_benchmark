#!/usr/bin/env python
import ConfigParser
import math
import multiprocessing
import operator
import os
import sys
import uuid
import zipfile

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

import cistargetx.common.progressmonitor as progressmonitor
import cistargetx.pwmrankings.orderstatistics as orderstatistics
from cistargetx.pwmrankings.createrankings import load_gene2regions


# Constants.
PROGRESS_FILE_PREFIX = "processed-motifs-"
PROGRESS_FILE_EXTENSION = ".lst"


def compute_rankings(cfg, pwms, species, outputfolder, log_output_fh, progress_output_fh):
    if cfg.has_option('regions', 'lut_file') and cfg.get('regions', 'lut_file'):
        gene2regions = load_gene2regions(cfg.get('regions', 'lut_file'))
    else:
        gene2regions = None

    os.nice(cfg.getint('multiprocessing', 'nice_increment'))
    pm = progressmonitor.ProgressMonitor(len(pwms), 12)

    for pwm, zipfilename in pwms:
        pm.start_iteration()
        log(cfg, "running orderstatistics for " + pwm + " " + pm.progress_str, log_output_fh)
        ranking_filenames = open_archive(zipfilename, species, outputfolder)
        filename = os.path.join(outputfolder, pwm + ".r")

        if len(ranking_filenames) != len(species):
            print >> sys.stderr, 'Not all wanted species "*.<species>.rr" files where found inside {0:s}.'.format(zipfilename)
            sys.exit(1)
        elif len(ranking_filenames) == 1:
            region2rank = dict()

            with open(ranking_filenames[0], 'r') as input_fh:
                for line in input_fh:
                    region, rank = line.rstrip().split('\t')
                    region2rank[region] = float(rank)

            # Mapping of regions to genes if necessary.
            if gene2regions:
                gene2rank = dict()
                for gene in gene2regions.keys():
                    gene2rank[gene] = min(region2rank[region] for region in gene2regions[gene])

                ranked_gene_or_region_ids = [gene for gene, score in sorted(gene2rank.items(),
                                                                            reverse=False,
                                                                            key=operator.itemgetter(1))]
            else:
                ranked_gene_or_region_ids = [region for region, score in sorted(region2rank.items(),
                                                                                reverse=False,
                                                                                key=operator.itemgetter(1))]
            with open(filename, 'w') as output_fh:
                for idx, gene_or_region_id in enumerate(ranked_gene_or_region_ids):
                    print >> output_fh, "{0:s}\t{1:d}".format(gene_or_region_id, idx + 1)
        else:
            with open(filename, 'w') as output_fh:
                for gene_or_region_id, rank in sorted(
                        orderstatistics.combined_rankings_iterator(cfg, ranking_filenames, gene2regions),
                        reverse=False,
                        key=operator.itemgetter(0)):
                    print >> output_fh, "{0:s}\t{1:d}".format(gene_or_region_id, rank)

        close_archive(ranking_filenames)

        print >> progress_output_fh, pwm
        progress_output_fh.flush()

        pm.end_iteration()

    if log_output_fh:
        log_output_fh.close()

    progress_output_fh.close()


def open_archive(zipfilename, species, outputfolder):
    filenames = []

    with zipfile.ZipFile(zipfilename, mode='r', allowZip64=True) as zip_input_fh:
        for filename in zip_input_fh.namelist():
            if not required_species(filename, species):
                continue

            zip_input_fh.extract(filename, path=outputfolder)
            fullfilename = os.path.join(outputfolder, filename)
            filenames.append(fullfilename)

    return filenames


def required_species(filename, species):
    for s in species:
        if filename.endswith('.' + s + '.rr'):
            return True

    return False


def close_archive(filenames):
    for filename in filenames:
        os.remove(filename)


def parse_arguments():
    if len(sys.argv) not in (3, 4):
        print >> sys.stderr, "Wrong number of input arguments."
        print >> sys.stderr, "Usage: python {0:s} <inifile> <outputfolder>".format(sys.argv[0])
        print >> sys.stderr, "Usage: python {0:s} <inifile> <species> <outputfolder>".format(sys.argv[0])
        sys.exit(2)

    inifile = sys.argv[1]

    if not os.path.isfile(inifile):
        print >> sys.stderr, "Ini file doesn't exist."
        sys.exit(1)

    cfg = ConfigParser.RawConfigParser()
    cfg.read(inifile)

    if len(sys.argv) == 3:
        species = cfg.get('recombination', 'include_species').split(';')
        outputdir = sys.argv[2]
    else:
        species = sys.argv[2].split(';')
        outputdir = sys.argv[3]

    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    return cfg, species, outputdir


def get_processed_motifs(outputdir):
    result = set()

    for name in os.listdir(outputdir):
        if name.startswith(PROGRESS_FILE_PREFIX) and name.endswith(PROGRESS_FILE_EXTENSION):
            with open(os.path.join(outputdir, name), 'r') as input_fh:
                for line in input_fh:
                    result.add(line.rstrip())

    return result


def log(cfg, msg, log_output_fh=None):
    verbose = cfg.getboolean('log', 'verbose')
    multiprocessing_enabled = cfg.getint('multiprocessing', 'number_of_cpus') > 1
    msg = "process-ID = " + str(os.getpid()) + "; " + msg if multiprocessing_enabled else msg

    if verbose:
        print msg
        sys.stdout.flush()
    if log_output_fh:
        print >> log_output_fh, msg
        log_output_fh.flush()


def main():
    cfg, species, outputdir = parse_arguments()

    processed_pwms = get_processed_motifs(outputdir)
    processed_pwms_count = len(processed_pwms)
    archive_extension = cfg.get('rankings', 'archive_extension')
    ranking_folder = cfg.get('rankings', 'ranking_folder')

    def derive_pwm_name(filename):
        idx = filename.index(archive_extension)
        return filename[:idx - 1]

    pwms = [(derive_pwm_name(filename), os.path.join(ranking_folder, filename))
            for filename in os.listdir(ranking_folder)
            if filename.endswith(archive_extension)
            and derive_pwm_name(filename) not in processed_pwms]
    remaining_pwms_count = len(pwms)

    if processed_pwms_count > 0:
        log(cfg,
            "{0:d}# pwms are already processed. Calculating rankings for remaining {1:d}# pwms".format(
                processed_pwms_count,
                remaining_pwms_count)
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
        log_output_fh = open(os.path.join(outputdir,
                                          log_filename_prefix + "-" + str(idx) + "." + log_filename_extension),
                             'w')
        progress_output_fh = open(os.path.join(outputdir,
                                               PROGRESS_FILE_PREFIX + str(uuid.uuid1()) + PROGRESS_FILE_EXTENSION),
                                  'w')

        start_idx = idx * block_size
        end_idx = min((idx + 1) * block_size, len(pwms))
        pwm_partition = pwms[start_idx:end_idx]

        p = multiprocessing.Process(target=compute_rankings,
                                    args=(cfg, pwm_partition, species, outputdir, log_output_fh, progress_output_fh))
        p.start()
        processes.append(p)

    for idx in range(process_nr):
        processes[idx].join()


if __name__ == "__main__":
    main()
