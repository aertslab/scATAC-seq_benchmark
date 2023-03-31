import os
import sys

from cistargetx.common.rankingsdatabase import RankingsDatabase


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, "Usage: python {0:s} <input_db1> <input_db2> ... <input_dbN> <output_db>".format(command_name)


def derive_name(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def main():
    if len(sys.argv) < 4:
        print >> sys.stderr, "Wrong number of input arguments."
        display_usage(sys.stderr)
        sys.exit(2)

    input_filenames = sys.argv[1:-1]
    output_filename = sys.argv[-1]

    if os.path.exists(output_filename):
        print >> sys.stderr, "Database file already exists."
        display_usage(sys.stderr)
        sys.exit(2)

    databases = []

    for filename in input_filenames:
        if not os.path.exists(filename):
            print >> sys.stderr, "Database file doesn't exists."
            display_usage(sys.stderr)
            sys.exit(2)

        databases.append(RankingsDatabase.load(filename, derive_name(filename)))

    reduce(lambda x, y: x.add_with_optimized_dtype(y), databases).write(output_filename)


if __name__ == "__main__":
    main()
