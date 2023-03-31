import sys
import os
import re

from cistargetx.common.rankingsdatabase import RankingsDatabase


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, "Usage: python {0:s} <input_db> <output_db> <feature_pattern>".format(command_name)


def derive_name(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def main():
    if len(sys.argv) != 4:
        print >> sys.stderr, "Wrong number of input arguments."
        display_usage(sys.stderr)
        sys.exit(2)

    input_filename, output_filename, feature_pattern = sys.argv[1:4]

    if os.path.exists(output_filename):
        print >> sys.stderr, "Database file already exists."
        display_usage(sys.stderr)
        sys.exit(2)

    if not os.path.exists(input_filename):
        print >> sys.stderr, "Database file doesn't exists."
        display_usage(sys.stderr)
        sys.exit(2)

    db = RankingsDatabase.load(input_filename, derive_name(input_filename))
    include_features = [name for name in db.feature_ids if not re.match(feature_pattern, name)]

    # TODO: weird bug "could not convert BLOB to buffer"
    # db.reduce_feature_ids(include_features).write(output_filename)
    db.retain_features_with_optimized_dtype(include_features).write(output_filename)


if __name__ == "__main__":
    main()
