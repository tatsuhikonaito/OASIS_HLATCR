#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
from collections import Counter
from sklearn.preprocessing import QuantileTransformer

# Set directories and file paths
base_dir = ""
sample_list_file = ""
dict_info_files = {"TCR": "", "BCR": ""}
os.chdir(base_dir)

# Load sample list
samples = np.loadtxt(sample_list_file, dtype=str)

# Fixed parameters
cell = "BCR"
adjust = "norm"
weight_method = "raw_cdr3"
full = "full"

def flatten_concatenation(matrix):
    return [item for row in matrix for item in row]

df = pd.DataFrame(index=samples)
info_file = "scRNA/BCR_SHM/data/BCR_postQCed_info_noflip_SHM_heavy.light_20240902.txt.gz"
info = pd.read_csv(info_file, sep="\t", header=0, index_col=0)

if full == "full":
    info = info.dropna(subset=["IR_VJ_1_c_call", "IR_VDJ_1_c_call", "IR_VDJ_1_d_call"])

seen_celltypes = set()
for col_celltype in ["l1", "l2", "l3"]:
    celltypes = sorted(set(info[col_celltype]))
    
    for celltype in ["all"] + celltypes:
        if celltype in seen_celltypes:
            continue
        seen_celltypes.add(celltype)
        info_cell = info if celltype == "all" else info[info[col_celltype] == celltype]
        
        for sample in samples:
            info_cell_sample = info_cell[info_cell.ID == sample]
            
            if weight_method == "raw_cdr3":
                info_cell_sample = info_cell_sample.drop_duplicates(subset=["clone_id"])
            
            info_cell_sample_K = info_cell_sample[info_cell_sample.IR_VJ_1_j_call.str[2] == "K"]
            info_cell_sample_L = info_cell_sample[info_cell_sample.IR_VJ_1_j_call.str[2] == "L"]
            
            for chain, group in zip(["H", "K", "L"], [info_cell_sample, info_cell_sample_K, info_cell_sample_L]):
                df.loc[sample, f"mean_mu_freq.{cell}.{chain}.{weight_method}.{full}.{celltype}"] = group["mu_freq_L" if chain in ["K", "L"] else "mu_freq_H"].mean()
                df.loc[sample, f"std_mu_freq.{cell}.{chain}.{weight_method}.{full}.{celltype}"] = group["mu_freq_L" if chain in ["K", "L"] else "mu_freq_H"].std()
                df.loc[sample, f"cnt_mu_freq.{cell}.{chain}.{weight_method}.{full}.{celltype}"] = len(group["mu_freq_L" if chain in ["K", "L"] else "mu_freq_H"])
                df.loc[sample, f"cv_mu_freq.{cell}.{chain}.{weight_method}.{full}.{celltype}"] = df.loc[sample, f"std_mu_freq.{cell}.{chain}.{weight_method}.{full}.{celltype}"] / df.loc[sample, f"mean_mu_freq.{cell}.{chain}.{weight_method}.{full}.{celltype}"]


# For Plink
qt = QuantileTransformer(random_state=0, output_distribution="normal")
df.iloc[:, :] = qt.fit_transform(df)

df.insert(0, "FID", "0")
df.insert(1, "IID", df.index)

output_filename = f"phenotypes/mu_freq.{cell}.ranknorm.pheno.gz"
df.to_csv(output_filename, sep="\t", header=True, index=False)
