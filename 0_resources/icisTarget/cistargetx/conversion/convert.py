import getopt
import os
import re
import sys


def load_filename(filename):
    with open(filename, 'r') as input_fh:
        return input_fh.read()


# TODO: support for GMT input and output.
# TODO: support for human region-based i-cisTarget and its nearest gene association between genes and regions.
def display_usage(ouput_fh=sys.stdout, cmd=sys.argv[0]):
    print >> ouput_fh, """Usage: python {0:s} [-f <overlap_fraction>] [-i <fbgn|cg|symbol|bed>] [-F gmt|legacy ] [-d <"up5kb+fullTx"|"up5kb+5utr+intron1"|"up5kb+introns"|"UCSC-up5kb+5utr+intron1">] <filename1> ... <filenameN>""".format(
        cmd)
    print >> ouput_fh
    print >> ouput_fh, "Optional parameters: -f overlap fraction (default = 0.40)"
    print >> ouput_fh, """                     -i input type (default = "bed")"""
    print >> ouput_fh, """                     -F gene signature output format (default = "legacy")"""
    print >> ouput_fh, """                     -d delineation (default="up5kb+5utr+intron1")"""


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:i:d:F:")
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        display_usage(sys.stderr)
        sys.exit(2)

    if len(args) < 1:
        print >> sys.stderr, "Wrong number of input arguments."
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    if not os.path.exists(args[0]):
        print >> sys.stderr, "'{0:s}' does not exist.".format(args[0])
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    overlap_fraction = 0.40
    input_type = 'bed'
    delineation = 'up5kb+5utr+intron1'
    output_format = 'legacy'

    for o, a in opts:
        if o == "-h":
            display_usage();
            sys.exit()
        elif o == "-f":
            try:
                overlap_fraction = float(a)
            except ValueError:
                print >> sys.stderr, "Overlap fraction must be a number between 0 and 1."
                print >> sys.stderr
                display_usage(sys.stderr)
                sys.exit(2)
        elif o == "-i":
            input_type = a
        elif o == "-d":
            delineation = a
        elif o == "-F":
            output_format = a
        else:
            assert False, "Unhandled option"

    filenames = args

    if overlap_fraction < 0.0 or overlap_fraction > 1.0:
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    if input_type not in ('fbgn', 'symbol', 'bed', 'cg'):
        print >> sys.stderr, "Wrong input type."
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    if delineation not in ('up5kb+fullTx', 'up5kb+5utr+intron1', 'up5kb+introns', 'UCSC-up5kb+5utr+intron1'):
        print >> sys.stderr, "Wrong delineation type."
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    if output_format not in ('gmt', 'legacy'):
        print >> sys.stderr, "Wrong output type."
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    if output_format == 'legacy' and len(filenames) > 1:
        print >> sys.stderr, "Only one input file possible with legacy output type."
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    # Import statement delayed until all arguments are parsed because this import
    # statement triggers the loading of many BED files into memory.
    from tools import bed2regions, gene_ids2regions

    if input_type == 'bed':
        for filename in filenames:
            query_data = load_filename(filename)
            region_ids = bed2regions(query_data, overlap_fraction)
            name = os.path.splitext(os.path.basename(filename))[0]

            if output_format == 'legacy':
                print >> sys.stdout, "\n".join(region_ids)
            else:
                print >> sys.stdout, name + "\t\t" + ",".join(region_ids)
    else:
        for filename in filenames:
            query_data = load_filename(filename)
            region_ids = gene_ids2regions(re.split('\n', query_data), gene_id_type=input_type, delineation=delineation,
                                          fraction=overlap_fraction)
            name = os.path.splitext(os.path.basename(filename))[0]

            if output_format == 'legacy':
                print >> sys.stdout, "\n".join(region_ids)
            else:
                print >> sys.stdout, name + "\t\t" + ",".join(region_ids)


if __name__ == "__main__":
    main()
