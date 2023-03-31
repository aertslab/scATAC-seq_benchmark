import os
import sys

from cistargetx.common.rankingsdatabase import RankingsDatabase


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, '\nUsage:   python {0:s} <input_db> <output_db> <gene_or_region_ids_filename>\n'.format(command_name)
    print >> output_fh, 'Purpose: Make rankings database from an existing one, but only retain rankings'
    print >> output_fh, '         for gene or region IDs in the gene_or_region_ids_filename.\n'


def derive_name(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def main():
    if len(sys.argv) != 4:
        print >> sys.stderr, "ERROR: Wrong number of input arguments."
        display_usage(sys.stderr)
        sys.exit(2)

    input_filename, output_filename, gene_or_region_ids_filename = sys.argv[1:4]

    if os.path.exists(output_filename):
        print >> sys.stderr, 'ERROR: Output database file "{0:s}" already exists.'.format(output_filename)
        display_usage(sys.stderr)
        sys.exit(2)

    if not os.path.exists(input_filename):
        print >> sys.stderr, 'ERROR: Input database file "{0:s}" does not exists.'.format(input_filename)
        display_usage(sys.stderr)
        sys.exit(2)

    if not os.path.exists(gene_or_region_ids_filename):
        print >> sys.stderr, 'ERROR: Input filename "{0:s}" with list of gene or region IDs does not exists.'.format(
            gene_or_region_ids_filename)
        display_usage(sys.stderr)
        sys.exit(2)

    include_gene_or_region_ids = set()

    with open(gene_or_region_ids_filename, 'r') as fh:
        for line in fh:
            include_gene_or_region_ids.add(line.rstrip('\r\n'))

    db = RankingsDatabase.load(input_filename, derive_name(input_filename), gene_ids=include_gene_or_region_ids)

    db.reduce_gene_ids_with_optimized_dtype(include_gene_or_region_ids).write(output_filename)


if __name__ == "__main__":
    main()
