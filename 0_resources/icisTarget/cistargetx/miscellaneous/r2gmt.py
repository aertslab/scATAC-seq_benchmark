#!/usr/bin/env python
import getopt
import operator
import os
import re
import sys


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, "Usage: python {0:s} [-f 2col|1col] <rankings_filename> <increment> <max_size>".format(
        command_name)
    print >> output_fh, "If 'stdin' is used as name for the rankingsfile, standard input is used."


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:h")
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        display_usage(sys.stderr)
        sys.exit(2)

    input_format = "2col"

    for o, a in opts:
        if o == "-f":
            if a.lower() not in ('2col', '1col'):
                print >> sys.stderr, "'{0:s}' is not a valid input format.".format(a)
                display_usage(sys.stderr)
                sys.exit(1)

            input_format = a.lower()
        elif o == "-h":
            display_usage(); sys.exit()
        else:
            assert False, "Unhandled option"

    if len(args) != 3:
        print >> sys.stderr, "Wrong number of input arguments."
        display_usage(sys.stderr)
        sys.exit(2)

    rankings_filename = args[0]
    inc, max_n = map(int, args[1:])

    with open(rankings_filename, 'r') if rankings_filename != 'stdin' else sys.stdin as input_fh:
        if input_format == "2col":
            id2r = dict([re.split('[ \t]+', line.rstrip()) for line in input_fh])
            ranked_ids = map(operator.itemgetter(0), sorted(id2r.iteritems(), key=operator.itemgetter(1)))
        else:
            ranked_ids = [line.rstrip() for line in input_fh]

    start_idx = 0
    end_idx = inc

    while (end_idx - start_idx) <= max_n:
        name = "{0:d}-{1:d}".format(start_idx + 1, end_idx)
        print >> sys.stdout, name + "\t\t" + ",".join(ranked_ids[start_idx:end_idx])
        end_idx += inc


if __name__ == "__main__":
    main()
