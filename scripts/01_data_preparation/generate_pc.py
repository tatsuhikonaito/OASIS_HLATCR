#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
from scipy.stats import zscore
from sklearn.decomposition import PCA
from sklearn.preprocessing import QuantileTransformer

# Set directory and file names
base_dir = ""
sample_list = ""
dict_info_files = {"TCR": "", "BCR": ""}
os.chdir(base_dir)

# Define constants
tcrs = ["TRA", "TRB"]
bcrs = ["IGH", "IGK", "IGL"]
cells = ["TCR", "BCR"]
genes = ["V", "D", "J"]
ces = ["ce", "raw_cdr3"]
chains_dict = {"TCR": ["A", "B"], "BCR": ["H", "K", "L"]}

metadata_trust = pd.read_csv("data/metadata.scrna.234.txt", sep="\t").set_index("sample.id")
samples = np.loadtxt(sample_list, dtype=str)
metadata_trust = metadata_trust.loc[samples]

# Fixed parameters
adjust = "norm"
weight_method = "raw_cdr3"
full = "full"
filtered = "unfiltered"

for cell in cells:
    chains = chains_dict[cell]
    info = pd.read_csv(dict_info_files[cell], sep="\t", index_col=0)
    celltypes = sorted(set(info["l1"]).union(info["l2"], info["l3"]))

    for gene in genes:
        for celltype in celltypes:
            prefix = f"{cell}.{gene}.qc.functional.{weight_method}.{full}.{celltype}.{adjust}"
            df_template = pd.DataFrame(index=metadata_trust.index, columns=["IID"])
            df_template["FID"] = "0"
            df_template["IID"] = df_template.index
            
            try:
                usage = pd.read_csv(f"usage.88.146/{prefix}.bed.gz", sep="\t")
                usage = usage.loc[:, list(usage.columns[:4]) + list(metadata_trust.index)]
            except FileNotFoundError:
                continue
                                
            usage = usage[samples].T
            
            for chain in chains:
                usage_chain = usage.loc[:, usage.columns.str[2] == chain]
                if (usage_chain.sum(axis=1) == 0).sum() > 10 or usage_chain.shape[1] < 3:
                    continue
                
                n_components = 2 if cell == "BCR" and gene == "J" else 5
                usage_chain = zscore(usage_chain)
                pcs = pd.DataFrame(PCA(n_components=n_components).fit_transform(usage_chain),
                                   index=usage_chain.index,
                                   columns=[f"pc{i}" for i in range(1, n_components + 1)])
                
                # For Plink
                df = df_template.copy()
                df = df.join(pcs)
                
                qt = QuantileTransformer(random_state=0, output_distribution='normal')
                df[pcs.columns] = qt.fit_transform(df[pcs.columns])
                output_filename = f"phenotypes/pc.{prefix}.{filtered}.{chain}.ranknorm.pheno.gz"
                df.to_csv(output_filename, sep="\t", index=False)
