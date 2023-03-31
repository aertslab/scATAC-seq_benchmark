import numpy
import operator
import os
import sys
import zipfile

from cistargetx.common.rankingsdatabase import RankingsDatabase
from cistargetx.common.utils import derive_dtype


RANKING_EXTENSION = ".r"


def _get_motif_name(filename):
    return filename[:filename.rindex(RANKING_EXTENSION)]


def _ranking_iterator(input_fh):
    for line in input_fh.readlines():
        columns = line.rstrip().split('\t')

        yield columns[0], int(float(columns[1])) - 1


def write_database(database_filename, motif_ids, get_file_handle):
    gene_or_region_ids = sorted(gene_or_region_id
                                for gene_or_region_id, rank in _ranking_iterator(get_file_handle(motif_ids[0])))
    gene_or_region_ids_count = len(gene_or_region_ids)

    print >> sys.stderr, "Warning: input must be 1-based encoding of rank."
    print >> sys.stderr, "Number of motifs = {0:d}".format(len(motif_ids))
    print >> sys.stderr, "Number of genes/regions = {0:d}".format(gene_or_region_ids_count)

    if gene_or_region_ids_count <= 2 ** 15:
        print >> sys.stderr, "=> int16, 0-based encoding of rank"
    else:
        print >> sys.stderr, "=> int32, 0-based encoding of rank"

    dtype = derive_dtype(gene_or_region_ids_count)
    # Fill rankings array with zeros, but do not use numpy.zeros as it does not preallocate the whole array.
    # Create the rankings array in Fortran order as we fill the array column by column.
    rankings = numpy.full(shape=(gene_or_region_ids_count, len(motif_ids)), fill_value=0, dtype=dtype, order='F')
    col_idx = 0

    for motif_id in motif_ids:
        motif_id_fh = get_file_handle(motif_id)
        current_gene_or_region_id_count = len(sorted(gene_or_region_id
                                                     for gene_or_region_id, rank in _ranking_iterator(motif_id_fh)))
        motif_id_fh.close()

        if gene_or_region_ids_count != current_gene_or_region_id_count:
            print >> sys.stderr, "Error: Ranking file '{0:s}' doesn\'t have the same number of genes/regions as ranking file '{1:s}'.".format(
                motif_id + '.r', motif_ids[0] + '.r')
            sys.exit(2)

        motif_id_fh = get_file_handle(motif_id)
        rankings[:, col_idx] = map(operator.itemgetter(1),
                                   sorted(_ranking_iterator(motif_id_fh), key=operator.itemgetter(0)))
        motif_id_fh.close()
        col_idx += 1

    RankingsDatabase.create(os.path.basename(database_filename),
                            motif_ids,
                            gene_or_region_ids,
                            rankings).write(database_filename)


def create_motif_id2filename_lut_for_zip(zip_input_fh):
    return dict((_get_motif_name(filename), filename)
                for filename in zip_input_fh.namelist() if filename.endswith(RANKING_EXTENSION))


def create_motif_id2filename_lut_for_folder(rankings_folder):
    return dict((_get_motif_name(filename), os.path.join(rankings_folder, filename))
                for filename in os.listdir(rankings_folder) if filename.endswith(RANKING_EXTENSION))


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, "Usage: python {0:s} <rankings_folder>|<rankings_archive> <databasefilename>".format(
        command_name)


def main():
    if len(sys.argv) != 3:
        print >> sys.stderr, "Wrong number of input arguments."
        display_usage(sys.stderr)
        sys.exit(2)

    if os.path.isdir(sys.argv[1]):
        rankings_folder = sys.argv[1]
        motif_id2filename = create_motif_id2filename_lut_for_folder(rankings_folder)
        zip_input_fh = None

        def get_file_handle(motif_id):
            return open(motif_id2filename[motif_id])
    else:
        rankings_archive = sys.argv[1]
        zip_input_fh = zipfile.ZipFile(rankings_archive, mode='r', allowZip64=True)
        motif_id2filename = create_motif_id2filename_lut_for_zip(zip_input_fh)

        def get_file_handle(motif_id):
            return zip_input_fh.open(motif_id2filename[motif_id], 'r')

    motif_ids = list(motif_id2filename.keys())
    motif_ids.sort()
    motif_count = len(motif_ids)

    if not motif_count:
        print >> sys.stderr, "No motif rankings present in supplied folder."
        display_usage(sys.stderr)
        sys.exit(2)

    database_filename = sys.argv[2]

    if os.path.exists(database_filename):
        print >> sys.stderr, "Database file already exists."
        display_usage(sys.stderr)
        sys.exit(2)

    write_database(database_filename, motif_ids, get_file_handle)

    if zip_input_fh:
        zip_input_fh.close()


if __name__ == "__main__":
    main()
