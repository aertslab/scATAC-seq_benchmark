#!/usr/bin/env python
import re
import sqlite3
import sys


def change_motif_id(motif_id):
    old2new_source_name = {'transfac10.1': 'transfac_pro', 'jaspar4': 'jaspar'}
    old_source_name, motif_name = re.split('-', motif_id.rstrip())[:2]

    if old_source_name in ("flyfactorsurvey", "selexconsensus"):
        return motif_id
    else:
        new_source_name = (old2new_source_name[old_source_name]
                           if old_source_name in old2new_source_name
                           else old_source_name)

        return "{0:s}-{1:s}".format(new_source_name, motif_name)


def main():
    if len(sys.argv) != 2:
        print >> sys.stderr, "Wrong number of input arguments."
        sys.exit(2)

    database_filename = sys.argv[1]

    with sqlite3.connect(database_filename) as connection:
        cursor = connection.cursor()
        cursor.execute('SELECT motifName, idx from motifs;')
        motif_ids = [row for row in cursor]

        print >> sys.stderr, "Nr of unique IDs (before correction) = {0:d}".format(len(motif_ids))

        for motif_id, idx in motif_ids:
            new_motif_id = change_motif_id(motif_id)
            print >> sys.stderr, "{0:d}:{1:s} => {2:s}".format(idx, motif_id, new_motif_id)
            cursor.execute('UPDATE motifs SET motifName = ? WHERE idx = ?;', (new_motif_id, idx))

        connection.commit()
        cursor.execute('SELECT motifName from motifs;')
        new_motif_ids = set(row[0] for row in cursor)
        print >> sys.stderr, "Nr of unique IDs (after correction) = {0:d}".format(len(new_motif_ids))


if __name__ == "__main__":
    main()
