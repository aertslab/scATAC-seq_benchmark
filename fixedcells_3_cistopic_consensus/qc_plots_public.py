import kde
import pycisTopic
import glob
import os
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import multiprocess as mp
from multiprocess import Pool
import pprint as pp

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


def threshold_otsu(array, nbins=100, min_value=100):
    """
    Apply Otsu threshold on topic-region distributions [Otsu, 1979].
    Parameters
    ---------
    array: `class::np.array`
            Array containing the region values for the topic to be binarized.
    nbins: int
            Number of bins to use in the binarization histogram
    Return
    ---------
    float
            Binarization threshold.
    Reference
    ---------
    Otsu, N., 1979. A threshold selection method from gray-level histograms. IEEE transactions on systems, man, and
    cybernetics, 9(1), pp.62-66.
    """
    array = array[(array >= min_value)]
    hist, bin_centers = histogram(array, nbins)
    hist = hist.astype(float)
    # Class probabilities for all possible thresholds
    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]
    # Class means for all possible thresholds
    mean1 = np.cumsum(hist * bin_centers) / weight1
    mean2 = (np.cumsum((hist * bin_centers)[::-1]) / weight2[::-1])[::-1]
    # Clip ends to align class 1 and class 2 variables:
    # The last value of ``weight1``/``mean1`` should pair with zero values in
    # ``weight2``/``mean2``, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2
    idx = np.argmax(variance12)
    threshold = bin_centers[:-1][idx]
    return threshold


# def calc_kde(xy):
#     return gaussian_kde(xy)(xy)


def plot_frag_qc(
    x,
    y,
    ax,
    x_thr_min=None,
    x_thr_max=None,
    y_thr_min=None,
    y_thr_max=None,
    ylab=None,
    xlab="Number of (unique) fragments",
    cmap="viridis",
    density_overlay=False,
    s=10,
    marker="+",
    c="#343434",
    xlim=None,
    ylim=None,
    **kwargs
):
    from scipy.stats import gaussian_kde

    assert all(x.index == y.index)
    barcodes = x.index.values

    if density_overlay:
        cores = 8

        x_log = np.log(x)

        # Split input array for KDE [log(x), y] array in
        # equaly spaced parts (start_offset + n * nbr_cores).
        kde_parts = [np.vstack([x_log[i::cores], y[i::cores]]) for i in range(cores)]

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
        x, y, z, barcodes = x_ordered[idx], y_ordered[idx], z[idx], barcodes[idx]
    else:
        z = c

    sp = ax.scatter(x, y, c=z, s=s, edgecolors=None, marker=marker, cmap=cmap, **kwargs)
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])

    # Start with keeping all barcodes.
    barcodes_to_keep = np.full(x.shape[0], True)

    # Filter barcodes out if needed based on thresholds:
    if x_thr_min is not None:
        ax.axvline(x=x_thr_min, color="r", linestyle="--")

    if x_thr_max is not None:
        ax.axvline(x=x_thr_max, color="r", linestyle="--")

    if y_thr_min is not None:
        ax.axhline(y=y_thr_min, color="r", linestyle="--")

    if y_thr_max is not None:
        ax.axhline(y=y_thr_max, color="r", linestyle="--")

    ax.set_xscale("log")
    ax.set_xmargin(0.01)
    ax.set_ymargin(0.01)
    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)

def plot_qc(
    sample,
    min_dict,
    max_dict,
    metadata_bc_df,
    include_kde=False,
    detailed_title=True,
    s=4,
    min_x_val=100,
    min_y_val=1
):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150)

    # calculate thresholds using a double otsu strategy:
    # first we calculate an otsu threshold using a minimum of 10 fragments
    # then, this otsu threshold is used as the minimum for the second iteration

    x_arr = np.log10(metadata_bc_df["Unique_nr_frag"])
    x_threshold_log = threshold_otsu(x_arr, nbins=5000, min_value=np.log10(min_x_val))
    x_threshold = 10**x_threshold_log

    y_arr = metadata_bc_df["TSS_enrichment"]
    y_threshold = threshold_otsu(y_arr, nbins=5000, min_value=min_y_val)

    # calculate cells passing filter
    metadata_bc_df_passing_filters = metadata_bc_df.loc[
        (metadata_bc_df.Unique_nr_frag > x_threshold)
        & (metadata_bc_df.TSS_enrichment > y_threshold)
    ]
    bc_passing_filters = metadata_bc_df_passing_filters.index

    # plot everything
    plot_frag_qc(
        x=metadata_bc_df["Unique_nr_frag"],
        y=metadata_bc_df["TSS_enrichment"],
        ylab="TSS Enrichment",
        s=s,
        x_thr_min=x_threshold,
        y_thr_min=y_threshold,
        xlim=[10, max_dict["Unique_nr_frag"]],
        ylim=[0, max_dict["TSS_enrichment"]],
        density_overlay=include_kde,
        ax=ax1,
    )
    plot_frag_qc(
        x=metadata_bc_df["Unique_nr_frag"],
        y=metadata_bc_df["FRIP"],
        x_thr_min=x_threshold,
        ylab="FRIP",
        s=s,
        xlim=[10, max_dict["Unique_nr_frag"]],
        ylim=[0, 1],
        density_overlay=include_kde,
        ax=ax2,
    )
    plot_frag_qc(
        x=metadata_bc_df["Unique_nr_frag"],
        y=metadata_bc_df["Dupl_rate"],
        x_thr_min=x_threshold,
        ylab="Duplicate rate per cell",
        s=s,
        xlim=[10, max_dict["Unique_nr_frag"]],
        ylim=[0, 1],
        density_overlay=include_kde,
        ax=ax3,
    )

    if detailed_title:
        med_nf = round(
            metadata_bc_df.loc[bc_passing_filters, "Unique_nr_frag"].median(), 2
        )
        med_tss = round(
            metadata_bc_df.loc[bc_passing_filters, "TSS_enrichment"].median(), 2
        )
        med_frip = round(metadata_bc_df.loc[bc_passing_filters, "FRIP"].median(), 2)
        title = f"{sample}: Kept {len(bc_passing_filters)} cells using Otsu filtering. Median Unique Fragments: {med_nf:.0f}. Median TSS Enrichment: {med_tss:.2f}. Median FRIP: {med_frip:.2f}\nUsed a minimum of {x_threshold:.2f} fragments and TSS enrichment of {y_threshold:.2f})"
    else:
        title = sample

    fig.suptitle(title, x=0.5, y=0.95, fontsize=10)
    return bc_passing_filters, fig


###################
# Run the loop here
###################

metadata_bc_pkl_list = sorted(glob.glob("cistopic_qc_out/*metadata_bc.pkl"))
metadata_bc_pkl_path_dict = {}
for metadata_bc_pkl_path in metadata_bc_pkl_list:
    sample = metadata_bc_pkl_path.split("/")[-1].split("__metadata_bc.pkl")[0]
    metadata_bc_pkl_path_dict[sample] = metadata_bc_pkl_path

# calculate global maxima to make nice equal plots:
metadata_bc_df_all = pd.DataFrame()
for sample in metadata_bc_pkl_path_dict.keys():
    with open(metadata_bc_pkl_path_dict[sample], "rb") as fh:
        metadata_bc_df = pickle.load(fh)
        metadata_bc_df_all = pd.concat([metadata_bc_df_all, metadata_bc_df])
        print(f"added {sample} metadata_bc_df")

max_dict = {}
min_dict = {}
for stat in metadata_bc_df_all.columns:
    max_dict[stat] = metadata_bc_df_all[stat].max()
    min_dict[stat] = metadata_bc_df_all[stat].min()

pp.pprint(min_dict)
pp.pprint(max_dict)

standard_min_x_val = 100
standard_min_y_val = 0.5
min_otsu_frags_dict = {}
min_otsu_tss_dict = {}
for metadata_bc_pkl_path in sorted(metadata_bc_pkl_path_dict.keys()):
    sample = metadata_bc_pkl_path.split("/")[-1].split("__metadata_bc.pkl")[0]
    tech = sample.split('_')[1]
    if tech == "ddseq":
        min_otsu_frags_dict[sample] = 200
    elif tech == "hydrop":
        min_otsu_frags_dict[sample] = 300
    else:
        min_otsu_frags_dict[sample] = standard_min_x_val

    if tech == "s3atac":
        min_otsu_tss_dict[sample] = 0
    else:
        min_otsu_tss_dict[sample] = standard_min_y_val

pp.pprint(min_otsu_frags_dict)

for sample in metadata_bc_pkl_path_dict.keys():
    pkl_path = f"selected_barcodes/{sample}_bc_passing_filters_otsu.pkl"
    if os.path.exists(pkl_path):
        print(f"{pkl_path} bc passing filters exists, skipping...")
    else:
        print(f"{pkl_path} bc passing filters does not exist yet, generating...")
        print(f"\tLoading {metadata_bc_pkl_path_dict[sample]}")
        with open(metadata_bc_pkl_path_dict[sample], "rb") as fh:
            metadata_bc_df = pickle.load(fh)

        print(f"\tFiltering cells and generating QC plots.")
        if not sample in min_otsu_frags_dict.keys():
            print(f"\t{sample} not in minimum dict! Using standard value of {standard_min_x_val}")
            min_x_val = standard_min_x_val
            min_y_val = standard_min_y_val
        else:
            min_x_val = min_otsu_frags_dict[sample]
            min_y_val = min_otsu_tss_dict[sample]

        bc_passing_filters, fig = plot_qc(sample=sample, min_dict=min_dict, max_dict=max_dict, metadata_bc_df=metadata_bc_df, include_kde=kde, min_x_val=min_x_val, min_y_val=min_y_val)

        plt.tight_layout()
        plt.savefig(f"plots_qc/{sample}_qc_otsu.png", dpi=300, facecolor="white")
        plt.show()
        plt.close()

        print(f"\tSaving...")
        with open(
            f"selected_barcodes/{sample}_bc_passing_filters_otsu.pkl", "wb"
        ) as fh:
            pickle.dump(bc_passing_filters, fh)
        fh.close()

        fh = open(f"selected_barcodes/{sample}_bc_passing_filters_otsu.txt", "w")
        for bc in list(bc_passing_filters):
            fh.write(bc + "\n")
        fh.close()

        metadata_bc_df.loc[bc_passing_filters].to_csv(f'selected_barcodes/{sample}_metadata_bc_df.tsv', sep='\t')
