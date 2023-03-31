import argparse
import kde
import os
import pandas as pd
import numpy as np
import multiprocess as mp
from multiprocess import Pool
import pprint as pp


def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description="Calculate Gaussian KDE in 2 dimensions.",)
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='tab-separated dataframe with 2 columns (x and y, include string header "x" and "y" for each column)')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Path to save KDE values')
    parser.add_argument('--ncores', '-n', type=str, required=True,
                        help='number of cores')

    return parser

def histogram(array, nbins=100):
    """
    Draw histogram from distribution and identify centers.
    Parameters
    ---------
    array: `class::np.array`
            Scores distribution
    nbins: int
            Number of bins to use in the histogram
    Return
    ---------
    float
            Histogram values and bin centers.
    """
    array = array.ravel().flatten()
    hist, bin_edges = np.histogram(array, bins=nbins, range=None)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    return hist, bin_centers

def main():
    """
    The main executable function
    """

    parser = make_argument_parser()
    args = parser.parse_args()
    in_path = args.input
    out_path = args.output
    cores = int(args.ncores)

    df = pd.read_csv(in_path, sep='\t', header=0)

    x = df["x"]
    y = df["y"]

    # Split input array for KDE [log(x), y] array in
    # equaly spaced parts (start_offset + n * nbr_cores).
    kde_parts = [np.vstack([x[i::cores], y[i::cores]]) for i in range(cores)]

    # Get nultiprocess context object to spawn processes.
    # mp_ctx = mp.get_context("spawn")

    # Calculate KDE in parallel.
    with Pool(processes=cores) as pool:
        results = pool.map(kde.calc_kde, kde_parts)

    z = np.concatenate(results)

    # now order x and y in the same way that z was ordered, otherwise random z value is assigned to barcode:

    x_ordered = np.concatenate([x[i::cores] for i in range(cores)])
    y_ordered = np.concatenate([y[i::cores] for i in range(cores)])

    idx = (
        z.argsort()
    )  # order based on z value so that highest value is plotted on top, and not hidden by lower values
    x, y, z = x_ordered[idx], y_ordered[idx], z[idx]

    pd.DataFrame({"x": x, "y": y, "kde": z}).to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
