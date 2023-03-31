#!/usr/bin/env python
import getopt
import math
import multiprocessing
import operator
import os
import sys

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

from cistargetx.cistargetdb.bunch import Bunch
from cistargetx.common.rankingsdatabase import RankingsDatabase


ORDERSTATISTICS_CMD = os.path.join(os.path.dirname(__file__), 'OrderStatistics.jar')
RANKINGS_FOLDER = os.path.expandvars('${DATADIR}/cistargetx/rankings/')


class InvalidFileFormatException(Exception):
    def __init__(self, value):
        Exception.__init__(self)
        self.value = value

    def __str__(self):
        return repr(self.value)


class GeneTree:
    @staticmethod
    def load_from_files(filenames):
        species2id2orthologs = dict()
        ref = None

        for filename in filenames:
            cur_ref, cur_dest = GeneTree._split_filename(filename)

            if not ref:
                ref = cur_ref
            elif ref != cur_ref:
                raise InvalidFileFormatException("Gene tree files are not compatible.")

            species2id2orthologs[cur_dest] = GeneTree._load_gene_tree(filename)

        return GeneTree(ref, species2id2orthologs)

    @staticmethod
    def _split_filename(filename):
        return os.path.splitext(os.path.basename(filename))[0].split('-')[:2]

    @staticmethod
    def _load_gene_tree(filename):
        id2orthologs = dict()

        with open(filename, 'r') as input_fh:
            for line in input_fh:
                columns = line.rstrip().split('\t')

                if len(columns) > 1 and len(columns[1].strip()) != 0:
                    id2orthologs.setdefault(columns[0], []).append(columns[1])

        return id2orthologs

    def __init__(self, ref_species, species2id2orthologs):
        self.ref_species = ref_species
        self.species2id2orthologs = species2id2orthologs

    def find_orthologs(self, species, gene_id):
        if gene_id in self.species2id2orthologs[species].keys():
            return self.species2id2orthologs[species][gene_id]
        else:
            return []


class RankingTable:
    @staticmethod
    def load_from_file(species, filename):
        name = os.path.splitext(os.path.basename(filename))[0]
        db = RankingsDatabase.load(filename, name)

        return RankingTable(species, db.feature_ids, db.gene_ids, db.rankings, db.total_gene_count, db.dtype)

    def __init__(self, species, motifs, gene_ids, rankings, totalgenecount, dtype):
        self.species = species
        self.motifs = motifs
        self.gene_ids = gene_ids
        self.rankings = rankings
        self.totalgenecount = totalgenecount
        self.dtype = dtype

    def get_rank_ratios(self, motif):
        col_idx = self.motifs.index(motif)
        rank_ratios = (self.rankings[:, col_idx] + 1.0) / self.totalgenecount

        return dict([(gene_id, rank_ratio) for gene_id, rank_ratio in zip(self.gene_ids, rank_ratios)])


def translate(ref_species, motif, species2tables, genetree):
    gene_ids = species2tables[ref_species].gene_ids
    destination_species = [other_species for other_species in species2tables.keys() if other_species != ref_species]
    gene_id2rrs = dict([(gene_id, [rr])
                        for gene_id, rr in species2tables[ref_species].get_rank_ratios(motif).iteritems()])

    for other_species in destination_species:
        for gene_id, rr in orthologous_rank_ratio_iterator(gene_ids, other_species, species2tables[other_species],
                                                           genetree, motif):
            gene_id2rrs.setdefault(gene_id, []).append(rr)

    return gene_id2rrs


def orthologous_rank_ratio_iterator(gene_ids, other_species, rankingtable, genetree, motif):
    totalgenecount = len(gene_ids)
    lost_gene_ids = set()

    def score_iterator(gene_ids, other_species, rankingtable, genetree, motif):
        orthologs_gene_id2rr = rankingtable.get_rank_ratios(motif)

        for gene_id in gene_ids:
            scores = [orthologs_gene_id2rr[orthologs_gene_id]
                      for orthologs_gene_id in genetree.find_orthologs(other_species, gene_id)
                      if orthologs_gene_id in orthologs_gene_id2rr.keys()]

            if len(scores) > 0:
                yield gene_id, min(scores)
            else:
                lost_gene_ids.add(gene_id)

    for idx, gene_id in enumerate(map(operator.itemgetter(0),
                                      sorted(score_iterator(gene_ids, other_species, rankingtable, genetree, motif),
                                             key=operator.itemgetter(1)))):
        yield gene_id, (idx + 1.0) / totalgenecount

    nan = float('nan')

    for gene_id in lost_gene_ids:
        yield gene_id, nan


def combined_rank_iterator(filename):
    exec_command = ['java', '-jar', ORDERSTATISTICS_CMD, filename]

    try:
        stdout, stderr = subprocess.Popen(exec_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    except OSError, msg:
        print >> sys.stderr, "orderstatistics: error during execution of '" + ' '.join(exec_command) + "': " + str(msg)
        sys.exit(1)

    if stderr and stderr != "":
        print >> sys.stderr, "orderstatistics: error during execution of '" + ' '.join(exec_command) + "': " + str(
            stderr)
        sys.exit(1)

    gene_id2rr = dict()

    for line in stdout.split('\n'):
        columns = line.split('\t')

        # Skip lines that don't contain a ranking.
        if len(columns) != 2:
            continue

        gene_id = columns[0]
        rr = float(columns[1])

        # CAVE: orderstatistics outputs only one ranking for each ID.
        gene_id2rr[gene_id] = rr

    ranked_gene_ids = [gene_id for gene_id, rr in sorted(gene_id2rr.items(), key=operator.itemgetter(1))]

    for idx, gene_id in enumerate(ranked_gene_ids):
        yield gene_id, idx + 1


def analyze_motifs(motifs, outputfolder, ref_species, species2tables, genetree):
    for idx, motif in enumerate(motifs):
        print >> sys.stderr, str(idx + 1) + ". " + motif + "(" + ref_species + ")"
        tblfilename = os.path.join(outputfolder, motif + ".tsv")

        with open(tblfilename, 'w') as output_fh:
            for gene_id, rrs in sorted(translate(ref_species, motif, species2tables, genetree).iteritems(),
                                       key=operator.itemgetter(0)):
                def to_str(rank):
                    return "{0:.16g}".format(rank)

                print >> output_fh, gene_id + "\t" + "\t".join(map(to_str, rrs))

        rfilename = os.path.join(outputfolder, motif + ".r")

        with open(rfilename, 'w') as output_fh:
            for gene_id, rank in combined_rank_iterator(tblfilename):
                print >> output_fh, gene_id + "\t" + str(rank)


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, "Usage: python {0:s} [-n <n_cpus>] [-d <primary_databases_folder>] <ref_species> <gene_trees_folder> <outputdir>".format(
        command_name)
    print >> output_fh
    print >> output_fh, "Optional parameters: -n <n_cpus>: number of cpus to use for calculations (default 12)"
    print >> output_fh, "                     -d <primary_databases_folder>: the folder containing the primary database (default ${DATADIR}/cistargetx/rankings/)"


def parse_arguments(args):
    try:
        opts, args = getopt.getopt(args, "hn:d:")
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        display_usage(sys.stderr)
        sys.exit(2)

    if len(args) != 3:
        print >> sys.stderr, "Wrong number of input arguments."
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    ref_species, gene_trees_foldername, output_foldername = args

    if not os.path.exists(gene_trees_foldername):
        print >> sys.stderr, "Invalid folder for gene trees: {0:s}.".format(gene_trees_foldername)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    n_cpus = 12
    primary_dbs_foldername = RANKINGS_FOLDER

    for o, a in opts:
        if o == "-h":
            display_usage();
            sys.exit()
        elif o == "-n":
            n_cpus = int(a)
        elif o == "-d":
            primary_dbs_foldername = a
        else:
            assert False, "Unhandled option"

    params = Bunch(n_cpus=n_cpus, primary_dbs_foldername=primary_dbs_foldername)

    return params, ref_species, gene_trees_foldername, output_foldername


def main():
    # Parsing input arguments.
    params, ref_species, gene_trees_foldername, output_foldername = parse_arguments(sys.argv[1:])

    if not os.path.exists(output_foldername):
        os.mkdir(output_foldername)

    # Loading data into memory.
    print >> sys.stderr, "Loading gene trees ..."
    genetree = GeneTree.load_from_files([os.path.join(gene_trees_foldername, filename)
                                         for filename in os.listdir(gene_trees_foldername)
                                         if filename.startswith(ref_species)])

    print >> sys.stderr, "Loading base rankings ..."
    species2tables = dict()

    for filename in os.listdir(params.primary_dbs_foldername):
        if filename.startswith('flybase') and filename.endswith('.db'):
            print >> sys.stderr, filename
            species = filename.split('-')[1]
            species2tables[species] = RankingTable.load_from_file(species,
                                                                  os.path.join(params.primary_dbs_foldername, filename))

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

    # Performing analysis.
    print >> sys.stderr, "Analyzing ..."
    motifs = species2tables[ref_species].motifs
    processes = []
    block_size = int(math.ceil(float(len(motifs)) / float(params.n_cpus)))

    for idx in range(params.n_cpus):
        start_idx, end_idx = idx * block_size, min((idx + 1) * block_size, len(motifs))
        motif_partition = motifs[start_idx:end_idx]

        p = multiprocessing.Process(target=analyze_motifs,
                                    args=(motif_partition, output_foldername, ref_species, species2tables, genetree))
        p.start()
        processes.append(p)

    for idx in range(params.n_cpus):
        processes[idx].join()


if __name__ == "__main__":
    main()
