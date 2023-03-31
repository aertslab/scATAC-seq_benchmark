#!/usr/bin/env python
import numpy
import os
import random
import sys

from databaseaccess import fetch_rankings, write_rankings


# TODO:


def main():
    # Parsing input arguments ...
    if len(sys.argv) != 4:
        print >> sys.stderr, "Wrong number of input arguments."
        print >> sys.stderr, "Usage: python {0:s} <n> <inputdbfilename> <outputdbfilename>".format(sys.argv[0])
        sys.exit(2)

    n = numpy.int32(sys.argv[1])

    input_database_filename, output_database_filename = sys.argv[2:]

    if os.path.exists(output_database_filename):
        print >> sys.stderr, "Database file already exists."
        sys.exit(2)

    # Loading data.
    features, gene_ids, rankings, totalgenecount, dtype = fetch_rankings(input_database_filename)
    featurecount, geneidcount = len(features), len(gene_ids)

    print >> sys.stderr, "Data loaded in memory ..."
    print >> sys.stderr, "Memory footprint rankings array = {0:d} bytes".format(rankings.itemsize * rankings.size)

    # Translate rankings matrix to geneid matrix.
    ranked_gene_ids = numpy.full(shape=(geneidcount, featurecount), fill_value=0, dtype=numpy.int_)
    rankingnrs = numpy.arange(geneidcount)

    for col_idx in range(featurecount):
        ranked_gene_ids[:, col_idx] = numpy.argsort(rankings[:, col_idx])

    print >> sys.stderr, "Memory representation changed ..."
    print >> sys.stderr, "Memory footprint gene_ids array = {0:d} bytes".format(
        ranked_gene_ids.itemsize * ranked_gene_ids.size)

    # Create random rankings.
    random_features = numpy.array(['random{0:d}'.format(idx) for idx in range(n)])
    # Fill rankings array with zeros, but do not use numpy.zeros as it does not preallocate the whole array.
    random_rankings = numpy.full(shape=(totalgenecount, n), fill_value=0, dtype=dtype)

    print >> sys.stderr, "Memory footprint random rankings array = {0:d} bytes".format(
        random_rankings.itemsize * random_rankings.size)
    for col_idx in range(n):
        # Fill random gene ID numbers array with zeros, but do not use numpy.zeros as it does not preallocate the whole array.
        random_geneid_nrs = numpy.full(shape=geneidcount, fill_value=0, dtype=numpy.int_)
        # Hashing as membership operation.
        previous_geneid_nrs = set()
        threshold = (1 * totalgenecount) / 3

        for row_idx in range(threshold):
            geneid_nr = -1

            while geneid_nr < 0 or geneid_nr in previous_geneid_nrs:
                # Draw uniform from ranked_gene_ids at given ranking position.
                geneid_nr = random.choice(ranked_gene_ids[row_idx, :])

            previous_geneid_nrs.add(geneid_nr)
            random_geneid_nrs[row_idx] = geneid_nr

        leftover = numpy.setdiff1d(numpy.arange(totalgenecount), random_geneid_nrs[:threshold])
        numpy.random.shuffle(leftover)
        random_geneid_nrs[threshold:] = leftover
        random_rankings[:, col_idx] = numpy.argsort(gene_ids[random_geneid_nrs])

        print >> sys.stderr, "\r{0:2d}".format(col_idx + 1),

    print >> sys.stderr, "Randomization finished ..."

    # Save random rankings.
    write_rankings(random_features, gene_ids, random_rankings, output_database_filename)

    print >> sys.stderr, "New database created and saved ..."


if __name__ == "__main__":
    main()
