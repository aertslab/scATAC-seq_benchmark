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
import kde

metadata_bc_pkl_path_dict = {
    x.split("/")[-1].split(f"__")[0]: x
    for x in sorted(
        glob.glob("/lustre1/project/stg_00090/scatac_benchmark/fixedcells_2_cistopic/cistopic_qc_out/*FIXEDCELLS*metadata*.pkl")
    )
}
print(metadata_bc_pkl_path_dict)

for sample in [
    "VIB_hydrop_11.FIXEDCELLS",
    "VIB_hydrop_12.FIXEDCELLS",
    "VIB_hydrop_21.FIXEDCELLS",
    "VIB_hydrop_22.FIXEDCELLS",
]:
    metadata_bc_pkl_path_dict.pop(sample)


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

metadata_bc_df_merged = pd.DataFrame()
metadata_bc_df_filtered_merged = pd.DataFrame()
df = pd.DataFrame(index=metadata_bc_pkl_path_dict.keys())
df["tech"] = [x.split("_")[1] for x in df.index]
techs = sorted(df["tech"].unique())

selected_bc_path_dict = {
    x.split("/")[-1].split(f"_bc")[0]: x
    for x in sorted(
        glob.glob("/lustre1/project/stg_00090/scatac_benchmark/fixedcells_2_cistopic/selected_barcodes/*FIXEDCELLS*.pkl")
    )
}
selected_bc_path_dict

for sample in [
    "VIB_hydrop_11.FIXEDCELLS",
    "VIB_hydrop_12.FIXEDCELLS",
    "VIB_hydrop_21.FIXEDCELLS",
    "VIB_hydrop_22.FIXEDCELLS",
]:
    selected_bc_path_dict.pop(sample)

tech_order = [
    "10xv2",
    "10xv1",
    "10xv11",
    "10xmultiome",
    "mtscatac",
    "ddseq",
    "s3atac",
    "hydrop",
]


for tech in techs:
    print(tech)
    samples = list(df[df["tech"] == tech].index)

    for sample in samples:
        metadata_path = metadata_bc_pkl_path_dict[sample]

        with open(metadata_path, "rb") as fh:
            metadata_bc_df = pickle.load(fh)

        metadata_bc_df["tech"] = sample.split("_")[1]

        with open(selected_bc_path_dict[sample], "rb") as fh:
            selected_bc = pickle.load(fh)

        selected_bc = [x.replace("FULL", "FIXEDCELLS") for x in selected_bc]

        metadata_bc_df_filtered = metadata_bc_df.loc[selected_bc]

        metadata_bc_df_merged = pd.concat([metadata_bc_df_merged, metadata_bc_df])
        metadata_bc_df_filtered_merged = pd.concat(
            [metadata_bc_df_filtered_merged, metadata_bc_df_filtered]
        )

metadata_bc_df_filtered_merged["Unique_nr_frag_in_regions"]

metadata_bc_df_filtered_merged["Log_unique_nr_frag_in_regions"] = np.log10(
    metadata_bc_df_filtered_merged["Unique_nr_frag_in_regions"]
)

metadata_bc_df_filtered_merged["tech_sample_id"] = (
    metadata_bc_df_filtered_merged["tech"]
    + "_"
    + metadata_bc_df_filtered_merged["sample_id"]
)

metadata_bc_df_filtered_merged["run"] = [
    x.split(".")[0][-1] for x in metadata_bc_df_filtered_merged["sample_id"]
]

metadata_bc_df_filtered_merged["replicate"] = [
    x.split(".")[0][-1] for x in metadata_bc_df_filtered_merged["sample_id"]
]

variables_list = [["Unique_nr_frag", "TSS_enrichment"], ["Unique_nr_frag", "FRIP"]]

n_rows = len(variables_list)
n_cols = len(metadata_bc_df_merged["tech"].unique())

# initiate axes
f, axes = plt.subplots(
    n_rows,
    n_cols,
    figsize=(n_cols * 3, n_rows * 3),
    # sharex="col",
    sharey="row",
    dpi=300,
)

tech_alias_dict = {
    "10xmultiome": "10x Multiome",
    "10xv1": "10x scATAC v1",
    "10xv11": "10x scATAC v1.1",
    "10xv2": "10x scATAC v2",
    "ddseq": "BioRad SureCell ATAC-seq",
    "hydrop": "HyDrop",
    "mtscatac": "mtscATAC-seq",
    "s3atac": "s3-ATAC",
}
var_alias_dict = {
    "Log_total_nr_frag": "Total Fragments",
    "Log_unique_nr_frag": "Total Fragments",
    "Total_nr_frag": "Total Fragments",
    "Unique_nr_frag": "Unique Fragments",
    "Dupl_nr_frag": "Duplicate Fragments",
    "Dupl_rate": "% Duplicate Fragments",
    "Total_nr_frag_in_regions": "Total Fragments in Regions",
    "Unique_nr_frag_in_regions": "Unique Fragments in Regions",
    "FRIP": "Fraction of Unique\nFragments in Peaks",
    "TSS_enrichment": "TSS Enrichment",
    "sample_id": "Sample",
    "tech": "Technology",
}

sample_order = [
    "STA_10xv11_1.FIXEDCELLS",
    "STA_10xv11_2.FIXEDCELLS",
    "TXG_10xv11_1.FIXEDCELLS",
    "VIB_10xv1_1.FIXEDCELLS",
    "VIB_10xv1_2.FIXEDCELLS",
    "CNA_10xv11_1.FIXEDCELLS",
    "CNA_10xv11_2.FIXEDCELLS",
    "CNA_10xv11_3.FIXEDCELLS",
    "CNA_10xv11_4.FIXEDCELLS",
    "CNA_10xv11_5.FIXEDCELLS",
    "VIB_10xv2_1.FIXEDCELLS",
    "VIB_10xv2_2.FIXEDCELLS",
    "CNA_10xv2_1.FIXEDCELLS",
    "CNA_10xv2_2.FIXEDCELLS",
    "TXG_10xv2_1.FIXEDCELLS",
    "TXG_10xv2_2.FIXEDCELLS",
    "HAR_ddseq_1.FIXEDCELLS",
    "HAR_ddseq_2.FIXEDCELLS",
    "UCS_ddseq_1.FIXEDCELLS",
    "UCS_ddseq_2.FIXEDCELLS",
    "BIO_ddseq_1.FIXEDCELLS",
    "BIO_ddseq_2.FIXEDCELLS",
    "BIO_ddseq_3.FIXEDCELLS",
    "BIO_ddseq_4.FIXEDCELLS",
    "SAN_10xmultiome_1.FIXEDCELLS",
    "SAN_10xmultiome_2.FIXEDCELLS",
    "VIB_10xmultiome_1.FIXEDCELLS",
    "VIB_10xmultiome_2.FIXEDCELLS",
    "CNA_10xmultiome_1.FIXEDCELLS",
    "CNA_10xmultiome_2.FIXEDCELLS",
    "BRO_mtscatac_1.FIXEDCELLS",
    "BRO_mtscatac_2.FIXEDCELLS",
    "CNA_mtscatac_1.FIXEDCELLS",
    "CNA_mtscatac_2.FIXEDCELLS",
    "MDC_mtscatac_1.FIXEDCELLS",
    "MDC_mtscatac_2.FIXEDCELLS",
    "OHS_s3atac_1.FIXEDCELLS",
    "OHS_s3atac_2.FIXEDCELLS",
    "VIB_hydrop_1.FIXEDCELLS",
    "VIB_hydrop_2.FIXEDCELLS",
    "EPF_hydrop_1.FIXEDCELLS",
    "EPF_hydrop_2.FIXEDCELLS",
    "EPF_hydrop_3.FIXEDCELLS",
    "EPF_hydrop_4.FIXEDCELLS",
    "CNA_hydrop_1.FIXEDCELLS",
    "CNA_hydrop_2.FIXEDCELLS",
    "CNA_hydrop_3.FIXEDCELLS",
]



cores = 36

include_kde = True
n_rows = 8
n_cols = 6
f, axes = plt.subplots(
    n_rows,
    n_cols,
    figsize=(n_cols * 3, n_rows * 3),
    sharex="col",
    sharey="row",
    dpi=300,
)


axes = axes.flatten()

include_kde = False
n_rows = 8
n_cols = 6
f, axes = plt.subplots(
    n_rows,
    n_cols,
    figsize=(n_cols * 3, n_rows * 3),
    sharex="col",
    sharey="row",
    dpi=300,
)

axes = axes.flatten()

for sample in sample_order:

    metadata_path = metadata_bc_pkl_path_dict[sample]

    with open(metadata_path, "rb") as fh:
        metadata_bc_df = pickle.load(fh)

    ## otsu here
    if not sample in min_otsu_frags_dict.keys():
        print(
            f"\t{sample} not in minimum dict! Using standard value of {standard_min_x_val}"
        )
        min_x_val = standard_min_x_val
        min_y_val = standard_min_y_val
    else:
        min_x_val = min_otsu_frags_dict[sample]
        min_y_val = min_otsu_tss_dict[sample]

    x_arr = np.log10(metadata_bc_df["Unique_nr_frag"])
    x_threshold_log = threshold_otsu(x_arr, nbins=5000, min_value=np.log10(min_x_val))
    x_threshold = 10**x_threshold_log

    y_arr = metadata_bc_df["TSS_enrichment"]
    y_threshold = threshold_otsu(y_arr, nbins=5000, min_value=min_y_val)

    ## plot
    ax = axes[sample_order.index(sample)]
    plot_qc(
        x=metadata_bc_df["Unique_nr_frag"],
        y=metadata_bc_df["TSS_enrichment"],
        ylab="TSS Enrichment",
        s=3,
        x_thr_min=x_threshold,
        y_thr_min=y_threshold,
        xlim=[10, max_dict["Unique_nr_frag"]],
        ylim=[0, max_dict["TSS_enrichment"]],
        density_overlay=include_kde,
        ax=ax,
    )

    ax.set_title(sample)

f.tight_layout()
plt.savefig("plts_general/LIBDS__frag_tss.png", facecolor="white", dpi=300)