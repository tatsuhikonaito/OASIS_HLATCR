#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
import numpy as np
from collections import Counter
from sklearn.preprocessing import QuantileTransformer

# Set directories and file paths
base_dir = ""
sample_list = ""
dict_info_files = {"TCR": "", "BCR": ""}
os.chdir(base_dir)

# Define constants
cells = ["TCR", "BCR"]

# Load sample list
samples = np.loadtxt(sample_list, dtype=str)

# Fixed parameters
adjust = "norm"
weight_method = "raw_cdr3"
full = "full"

# Process data for each cell type
for cell in cells:
    df = pd.DataFrame(index=samples)
    
    info = pd.read_csv(dict_info_files[cell], sep="\t", header=0, index_col=0)
    info["len_cdr3_1"] = info["IR_VJ_1_junction_aa"].str.len()
    info["len_cdr3_2"] = info["IR_VDJ_1_junction_aa"].str.len()
    
    if full == "full":
        info = info.dropna(subset=["IR_VJ_1_c_call", "IR_VDJ_1_c_call", "IR_VDJ_1_d_call"])
    
    seen_celltypes = set()    
    for col_celltype in ["l1", "l2", "l3"]:
        celltypes = sorted(set(info[col_celltype]))
        
        for celltype in celltypes:
            if celltype in seen_celltypes:
                continue
            seen_celltypes.add(celltype)
            
            info_cell = info[info[col_celltype] == celltype]
            
            for sample in samples:
                info_sample = info_cell[info_cell.ID == sample]
                
                if cell == "TCR":
                    if weight_method == "raw_cdr3":
                        info_sample = info_sample.drop_duplicates(subset=["clone_id"])
                    
                    for chain, col in zip(["A", "B"], ["len_cdr3_1", "len_cdr3_2"]):
                        df.loc[sample, f"mean_cdr3.{cell}.{chain}.{weight_method}.{full}.{celltype}"] = info_sample[col].mean()
                        df.loc[sample, f"std_cdr3.{cell}.{chain}.{weight_method}.{full}.{celltype}"] = info_sample[col].std()
                        df.loc[sample, f"cnt_cdr3.{cell}.{chain}.{weight_method}.{full}.{celltype}"] = len(info_sample[col])
                
                elif cell == "BCR":
                    if weight_method == "raw_cdr3":
                        info_sample = info_sample.drop_duplicates(subset=["clone_id"])
                    
                    info_sample_K = info_sample[info_sample.IR_VJ_1_j_call.str[2] == "K"]
                    info_sample_L = info_sample[info_sample.IR_VJ_1_j_call.str[2] == "L"]
                    
                    for chain, group in zip(["H", "K", "L", "KH", "LH"], [info_sample, info_sample_K, info_sample_L, info_sample_K, info_sample_L]):
                        col = "len_cdr3_2" if chain in ["H", "KH", "LH"] else "len_cdr3_1"
                        df.loc[sample, f"mean_cdr3.{cell}.{chain}.{weight_method}.{full}.{celltype}"] = group[col].mean()
                        df.loc[sample, f"std_cdr3.{cell}.{chain}.{weight_method}.{full}.{celltype}"] = group[col].std()
                        df.loc[sample, f"cnt_cdr3.{cell}.{chain}.{weight_method}.{full}.{celltype}"] = len(group[col])
                        df.loc[sample, f'cv_cdr3.{cell}.{chain}.{weight_method}.{full}.{celltype}'] = df.loc[sample, f"std_cdr3.{cell}.{chain}.{weight_method}.{full}.{celltype}"]/df.loc[sample, f"mean_cdr3.{cell}.{chain}.{weight_method}.{full}.{celltype}"]

    # For Plink
    qt = QuantileTransformer(random_state=0, output_distribution="normal")
    df.iloc[:, :] = qt.fit_transform(df)

    df.insert(0, "FID", "0")
    df.insert(1, "IID", df.index)
    output_filename = f"phenotypes/cdr3_len.{cell}.ranknorm.pheno.gz"
    df.to_csv(output_filename, sep="\t", header=True, index=False)
