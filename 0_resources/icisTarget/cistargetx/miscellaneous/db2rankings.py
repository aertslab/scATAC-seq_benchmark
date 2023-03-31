import os
import sys

from cistargetx.common.rankingsdatabase import RankingsDatabase


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, "Usage: python {0:s} <databasefilename> <feature_id1> ... <feature_idN>".format(command_name)


def main():
    if len(sys.argv) < 2:
        print >> sys.stderr, "Wrong number of input arguments."
        display_usage(sys.stderr)
        sys.exit(2)

    database_filename = sys.argv[1]

    if len(sys.argv) == 2:
        # Extract all feature IDs from the database.
        feature_ids = None
    else:
        feature_ids = sys.argv[2:]

    name = os.path.splitext(os.path.basename(database_filename))[0]
    db = RankingsDatabase.load(database_filename, name, feature_ids=feature_ids)

    if len(db.feature_ids) > 1:
        print >> sys.stdout, "#." + "\t".join(db.feature_ids)

    for idx, gene_id in enumerate(db.gene_ids):
        print >> sys.stdout, gene_id + "\t" + "\t".join(map(str, db.rankings[idx, :] + 1))


if __name__ == "__main__":
    main()
