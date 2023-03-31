#!/usr/bin/env python3

import gzip
import pickle
import argparse
import os
import numpy as np
import scrublet as scr
import scanpy as sc
from matplotlib import pyplot


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser(description='Template script')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input h5ad file containing the raw data'
)

parser.add_argument(
    "-o", "--output-prefix",
    dest="output_prefix",
    type=str,
    help='Output prefix of the files generated by this scripts.'
)

parser.add_argument(
    "-s", "--synthetic-doublet-umi-subsampling",
    type=float,
    dest="synthetic_doublet_umi_subsampling",
    default=1.0,
    help='Rate for sampling UMIs when creating synthetic doublets.'
)

parser.add_argument(
    "-m", "--min-counts",
    type=int,
    dest="min_counts",
    default=3,
    help='Used for gene filtering prior to PCA. Genes expressed at fewer than `min_counts` in fewer than `min_cells` (see below) are excluded.'
)

parser.add_argument(
    "-n", "--min-cells",
    type=int,
    dest="min_cells",
    default=3,
    help='Used for gene filtering prior to PCA. Genes expressed at fewer than `min_counts` (see above) in fewer than `min_cells` are excluded.'
)

parser.add_argument(
    "-v", "--min-gene-variability-pctl",
    type=float,
    dest="min_gene_variability_pctl",
    default=0.85,
    help='''
        Used for gene filtering prior to PCA. Keep the most highly variable genes (in the top min_gene_variability_pctl percentile),
        as measured by the v-statistic [Klein et al., Cell 2015].
        '''
)

parser.add_argument(
    "-l", "--log-transform",
    type=str2bool,
    action="store",
    dest="log_transform",
    default=False,
    help='''
        If True, log-transform the counts matrix (log10(1+TPM)).
        `sklearn.decomposition.TruncatedSVD` will be used for dimensionality reduction, unless `mean_center` is True.
        '''
)

parser.add_argument(
    "-c", "--mean-center",
    type=str2bool,
    action="store",
    dest="mean_center",
    default=True,
    help='''
        If True, center the data such that each gene has a mean of 0.  
        `sklearn.decomposition.PCA` will be used for dimensionality
        '''
)

parser.add_argument(
    "-w", "--normalize-variance",
    type=str2bool,
    action="store",
    dest="normalize_variance",
    default=True,
    help='''
If True, normalize the data such that each gene has a variance of 1.
`sklearn.decomposition.TruncatedSVD` will be used for dimensionality
reduction, unless `mean_center` is True.
        '''
)

parser.add_argument(
    "-p", "--n-prin-comps",
    type=int,
    dest="n_prin_comps",
    default=30,
    help='Number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph construction.'
)

parser.add_argument(
    "-r", "--threshold",
    type=float,
    dest="threshold",
    default=None,
    help="""
Doublet score threshold for calling a transcriptome a doublet. 
If `None`, this is set automatically by looking for the minimum between the two modes of the `doublet_scores_sim_` histogram. 
It is best practice to check the threshold visually using the `doublet_scores_sim_` histogram and/or based on co-localization of predicted doublets in a 2-D embedding."""
)

parser.add_argument(
    "-t", "--technology",
    type=str,
    dest="technology",
    choices=["10x"],
    help='Single-cell technology used.'
)

parser.add_argument(
    "-f", "--use-variable-features",
    type=str2bool,
    action="store",
    dest="use_variable_features",
    default=False,
    help='''
        If False, don't filter the feature space by the most variable featurse.
        '''
)

parser.add_argument(
    "-z", "--h5ad-with-variable-features-info",
    dest="h5ad_with_variable_features_info",
    type=argparse.FileType('r'),
    help='Input is a .h5ad containing the highly_variable slot. This should be coupled with the use of --use-variable-features True.'
)

args = parser.parse_args()


# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = args.output_prefix
SAMPLE_NAME = FILE_PATH_OUT_BASENAME.split(".")[0]

# I/O
# Expects h5ad file
try:
    adata_raw = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Wrong input format. Expects .h5ad files.")

################################################################################
# Processing...

if args.use_variable_features:
    print("Subsetting the variable features from the counts matrix...")
    if args.h5ad_with_variable_features_info is None:
        raise Exception("VSN ERROR: Expecting --h5ad-with-variable-features-info argument to be set since --use-variable-features argument is set to True.")

    FILE_PATH_H5AD_WITH_HVG_INFO = args.h5ad_with_variable_features_info
    adata_hvg = sc.read_h5ad(filename=FILE_PATH_H5AD_WITH_HVG_INFO.name)
    counts_matrix = adata_raw.X[:, np.array(adata_hvg.var['highly_variable'])]
else:
    counts_matrix = adata_raw.X

scrub = scr.Scrublet(counts_matrix)
adata_raw.obs['doublet_scores'], adata_raw.obs['predicted_doublets'] = scrub.scrub_doublets(
    synthetic_doublet_umi_subsampling=args.synthetic_doublet_umi_subsampling,
    use_approx_neighbors=True,
    distance_metric='euclidean',
    get_doublet_neighbor_parents=False,
    min_counts=args.min_counts,
    min_cells=args.min_cells,
    min_gene_variability_pctl=args.min_gene_variability_pctl,
    log_transform=args.log_transform,
    mean_center=args.mean_center,
    normalize_variance=args.normalize_variance,
    n_prin_comps=args.n_prin_comps,
    verbose=True
)


def save_histograms(out_basename, scrublet):
    fig, _ = scrublet.plot_histogram(fig_size=(20, 10))
    fig.suptitle(f"{out_basename} - Scrublet Histgrams", fontsize=16)
    fig.subplots_adjust(top=0.85)
    pyplot.savefig(f"{FILE_PATH_OUT_BASENAME}.ScrubletHistograms.png")


# Check if algorithm failed at identifying a doublet score threshold
if (adata_raw.obs["predicted_doublets"].isnull().values.any() or scrub.predicted_doublets_ is None) and args.threshold is None:
    # Set dummy threshold in order to save the plots
    scrub.threshold_ = 1
    # Save the histograms in the isolated process directory so that user can use it to define a doublet score threshold for the current sample.
    save_histograms(
        out_basename=FILE_PATH_OUT_BASENAME,
        scrublet=scrub
    )
    raise Exception(f"""
VSN ERROR: Scrublet failed to automatically identify a doublet score threshold for {SAMPLE_NAME}.
A manual doublet score threshold can be set using the --threshold (params.tools.scrublet.threshold) argument.
Consider to use sample-based parameter setting as described at https://vsn-pipelines.readthedocs.io/en/develop/features.html#multi-sample-parameters. E.g.:
params {{
    tools {{
        scrublet {{
            threshold = [
                {SAMPLE_NAME}: [your-custom-threshold-for-that-sample],
                ...
            ]
        }}
    }}
}}
In order to facilitate this process, the Scrublet histograms have been saved in the given process work directory.""")

# Call the doublets using custom threshold
if (adata_raw.obs["predicted_doublets"].isnull().values.any() or scrub.predicted_doublets_ is None) and args.threshold is not None:
    print(f"VSN MSG: Calling doublets using manual threshold: {args.threshold}")
    adata_raw.obs["predicted_doublets"] = scrub.call_doublets(threshold=args.threshold)


# Rename the columns
adata_raw.obs.rename(
    columns={
        "doublet_scores": "scrublet__doublet_scores",
        "predicted_doublets": "scrublet__predicted_doublets"
    },
    inplace=True
)

if args.technology == "10x":
    # Take doublet cells based on expected doublets based on number of cells (10x Chromium)
    cells_recovered = len(adata_raw)
    doublet_rate = 0.0008 * cells_recovered + 0.0527
    expected_doublets = np.int(doublet_rate / 100 * cells_recovered)
    doublet_cells = adata_raw.obs['scrublet__doublet_scores'].sort_values(
        ascending=False
    ).head(
        n=expected_doublets
    ).index
    adata_raw.obs['scrublet__predicted_doublets_based_on_10x_chromium_spec'] = False
    adata_raw.obs.loc[
        doublet_cells,
        'scrublet__predicted_doublets_based_on_10x_chromium_spec'
    ] = True
    final_adata = adata_raw.obs.loc[
        :,
        [
            'sample_id',
            'scrublet__doublet_scores',
            'scrublet__predicted_doublets',
            'scrublet__predicted_doublets_based_on_10x_chromium_spec'
        ]
    ]
else:
    raise Exception(f"Doublet detection with Scrublet for the given technolog {args.technology} is not implemented")

################################################################################

# I/O
final_adata.to_csv(
    path_or_buf=f"{FILE_PATH_OUT_BASENAME}.ScrubletDoubletTable.tsv",
    sep="\t",
    index=True,
    index_label='index',
    header=True
)
f = gzip.open(f"{FILE_PATH_OUT_BASENAME}.ScrubletObject.pklz", 'wb')
pickle.dump(scrub, f, protocol=4)
f.close()
