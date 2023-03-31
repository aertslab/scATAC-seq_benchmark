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


# def calc_kde(xy):
#     return gaussian_kde(xy)(xy)
metadata_bc_pkl_list = sorted(glob.glob("cistopic_qc_out/*metadata_bc.pkl"))
metadata_bc_pkl_path_dict = {}
for metadata_bc_pkl_path in metadata_bc_pkl_list:
    sample = metadata_bc_pkl_path.split("/")[-1].split("__metadata_bc.pkl")[0]
    metadata_bc_pkl_path_dict[sample] = metadata_bc_pkl_path


variables_list = [("Unique_nr_frag", "TSS_enrichment"), ("Unique_nr_frag", "FRIP")]

for sample in all_order:
    tech = sample.split("_")[1]
    pkl_path = f"selected_barcodes/{sample}_bc_passing_filters_otsu.pkl"
    print(f"{pkl_path} bc passing filters does not exist yet, generating...")
    print(f"\tLoading {metadata_bc_pkl_path_dict[sample]}")
    with open(metadata_bc_pkl_path_dict[sample], "rb") as fh:
        metadata_bc_df = pickle.load(fh)

    for x_var, y_var in variables_list:
        print(x_var)
        print(y_var)

        cores = 8
        x_log = np.log(metadata_bc_df[x_var])
        y = metadata_bc_df[y_var]

        kde_parts = [np.vstack([x_log[i::cores], y[i::cores]]) for i in range(cores)]

        with Pool(processes=cores) as pool:
            results = pool.map(kde.calc_kde, kde_parts)

        z = np.concatenate(results)

        x_ordered = np.concatenate([x[i::cores] for i in range(cores)])
        y_ordered = np.concatenate([y[i::cores] for i in range(cores)])

        idx = (
            z.argsort()
        )  # order based on z value so that highest value is plotted on top, and not hidden by lower values
        x, y, z, barcodes = x_ordered[idx], y_ordered[idx], z[idx], barcodes[idx]
        
        df = pd.DataFrame([x, y, z], columns=[x_var, y_var, f"{y_var}_density"])
        df.to_csv(f"intermediate_df/{sample}__{y_var}_density")
        

# metadata_bc_df.loc[bc_passing_filters].to_csv(f'selected_barcodes/{sample}_metadata_bc_df.tsv', sep='\t')