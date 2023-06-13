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

sns.set_context("notebook")
sns.set_style("darkgrid")

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
    "10xv1",
    "10xv11",
    "10xv2",
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

        # This line is new: creates a new dataframe that's metadata_bc_df minus the rows in selected_bc
        metadata_bc_df_minus_filtered = metadata_bc_df.drop(selected_bc, errors='ignore')

        metadata_bc_df_filtered = metadata_bc_df.loc[selected_bc]

        metadata_bc_df_merged = pd.concat([metadata_bc_df_merged, metadata_bc_df_minus_filtered]) # use the new df here
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



tech_alias_dict = {
    "10xmultiome": "10x MO",
    "10xv1": "10x v1",
    "10xv11": "10x v1.1",
    "10xv2": "10x v2",
    "ddseq": "ddSEQ SureCell",
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

tech_color_palette = {
    "10xv2": "#1b9e77",
    "10xv1": "#d95f02",
    "10xv11": "#7570b3",
    "10xmultiome": "#e7298a",
    "mtscatac": "#66a61e",
    "ddseq": "#e6ab02",
    "s3atac": "#a6761d",
    "hydrop": "#666666",
}
## load selected bc
selected_bc_path_dict = {
    x.split("/")[-1].split(f"_bc")[0]: x
    for x in sorted(
        glob.glob("/lustre1/project/stg_00090/scatac_benchmark/fixedcells_2_cistopic/selected_barcodes/*FIXEDCELLS*.pkl")
    )
}

print(selected_bc_path_dict)

for sample in [
    "VIB_hydrop_11.FIXEDCELLS",
    "VIB_hydrop_12.FIXEDCELLS",
    "VIB_hydrop_21.FIXEDCELLS",
    "VIB_hydrop_22.FIXEDCELLS",
]:
    selected_bc_path_dict.pop(sample)

selected_bc_path_dict

metadata_bc_df_merged = pd.DataFrame()
metadata_bc_df_filtered_merged = pd.DataFrame()

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

        # This line is new: creates a new dataframe that's metadata_bc_df minus the rows in selected_bc
        metadata_bc_df_minus_filtered = metadata_bc_df.drop(selected_bc, errors='ignore')

        metadata_bc_df_filtered = metadata_bc_df.loc[selected_bc]

        metadata_bc_df_merged = pd.concat([metadata_bc_df_merged, metadata_bc_df_minus_filtered]) # use the new df here
        metadata_bc_df_filtered_merged = pd.concat(
            [metadata_bc_df_filtered_merged, metadata_bc_df_filtered]
        )


cores = 32

# initiate axes
f, axes = plt.subplots(
    n_rows,
    n_cols,
    figsize=(n_cols * 3, n_rows * 3),
    # sharex="col",
    sharey="row",
    dpi=300,
)

# initiate axes
f, axes = plt.subplots(
    n_rows,
    n_cols,
    figsize=(n_cols * 3, n_rows * 3),
    # sharex="col",
    sharey="row",
    dpi=300,
)

print("plotting...")
naked = False
include_kde = True
scatterplots = []
for tech in tech_order:
    # subset df to tech
    df_tmp = metadata_bc_df_merged[metadata_bc_df_merged["tech"] == tech]
    print(tech)

    for variables in variables_list:
        x_var = variables[0]
        y_var = variables[1]
        ## calculate KDE

        x = df_tmp[x_var]
        x_log = np.log(df_tmp[x_var])
        y = df_tmp[y_var]
        barcodes = x.index.values

        if include_kde == True:
            # Split input array for KDE [log(x), y] array in
            # equaly spaced parts (start_offset + n * nbr_cores).
            kde_parts = [
                np.vstack([x_log[i::cores], y[i::cores]]) for i in range(cores)
            ]

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

            ## plot
            ax = axes[variables_list.index(variables), tech_order.index(tech)]
            scatter1 = sns.scatterplot(x=x, y=y, hue=z, palette="viridis", s=2, linewidth=0, ax=ax)

        else:
            z = None

            ## plot
            ax = axes[variables_list.index(variables), tech_order.index(tech)]
            scatter1 = sns.scatterplot(
                x=x,
                y=y,
                # hue=z,
                # palette="viridis",
                s=2,
                linewidth=0,
                ax=ax,
            )

        # plot cell barcodes on top
        df_tmp_filtered = metadata_bc_df_filtered_merged[
            metadata_bc_df_filtered_merged["tech"] == tech
        ]
        x = df_tmp_filtered[x_var]
        y = df_tmp_filtered[y_var]

        scatter2 = sns.scatterplot(
            x=x,
            y=y,
            color=tech_color_palette[tech],
            # palette="viridis",
            s=2,
            linewidth=0,
            ax=ax,
        )

        # ax.get_legend().remove()
        ax.set_xlim([10, max_dict[x_var]])
        ax.set_ylim([0, max_dict[y_var]])
        ax.set_xscale("log")

        # only set title on top row
        if variables == variables_list[0]:
            ax.set_title(tech_alias_dict[tech], fontsize=22)
        else:
            ax.set_title(None)

        # only set x label on bottom row
        if variables == variables_list[-1]:
            ax.set_xlabel(var_alias_dict[x_var], fontsize=14)
        else:
            ax.set_xlabel(None)

        # only set y label on left col
        if tech == tech_order[0]:
            ax.set_ylabel(var_alias_dict[y_var], fontsize=14)
        else:
            ax.set_ylabel(None)

        sns.despine(ax=ax)
        if naked == True:
            ax.set_xlim([10, max_dict[x_var]])
            ax.set_ylim([0, max_dict[y_var]])
            ax.set_xscale("log")

            ax.set_title(None)

            ax.set_xlabel(None)

            ax.set_ylabel(None)

            sns.despine(ax=ax, left=True, bottom=True)
            ax.axis("off")
            
        ax.get_legend().remove()
        scatterplots.extend([scatter1, scatter2])

plt.tight_layout()
f.tight_layout()

for scatter in scatterplots:
    scatter.set_visible(False)

plt.tight_layout()
plt.savefig("plts_general/FIXEDCELLS_distribution_bytech_kde_style_withaxes_v3_AXESONLY.svg", format='svg')


# Save PNG with only the scatterplot dots
for scatter in scatterplots:
    scatter.set_visible(True)

# Turn off the axis lines and labels
ax = plt.gca()
ax.axis('off')

plt.tight_layout()
plt.savefig("plts_general/FIXEDCELLS_distribution_bytech_kde_style_withaxes_v3_DOTSONLY.png", dpi=900)
