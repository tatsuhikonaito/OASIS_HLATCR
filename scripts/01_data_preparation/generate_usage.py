#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
from collections import Counter
from sklearn.preprocessing import QuantileTransformer

# Set directory and file names
base_dir = ""
sample_list_file = ""
dict_info_files = {"TCR": "", "BCR": ""}
os.chdir(base_dir)

# Define constants
cells = ["TCR", "BCR"]
dict_chains = {
    "TCR": ["TRA", "TRB"],
    "BCR": {"V": ["IGH", "IGK", "IGL"], "D": ["IGH"], "J": ["IGH", "IGK", "IGL"]}
}
genes = ["V", "D", "J"]
weight_methods = ["raw_cdr3", "weighted"]

# Load sample list
samples = np.loadtxt(sample_list_file, dtype=str)

# Fixed parameters
adjust = "norm"
full = "full"

def flatten_concatenation(matrix):
    return [item for row in matrix for item in row]

for cell in cells:
    info = pd.read_csv(dict_info_files[cell], sep="\t", header=0, index_col=0)

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
            
            for gene in genes:
                if gene == "D":
                    index = sorted(set(info[f"IR_VDJ_1_{gene.lower()}_call"].dropna()))
                else:
                    index = sorted(set(info[f"IR_VJ_1_{gene.lower()}_call"])) + \
                            sorted(set(info[f"IR_VDJ_1_{gene.lower()}_call"]))

                for weight_method in weight_methods:
                    usage = pd.DataFrame(0, columns=["#Chr", "start", "end", "ID"] + list(samples), index=index)

                    for sample in samples:
                        info_sample = info_cell[info_cell.ID == sample]

                        if weight_method == "raw_cdr3":
                            info_sample = info_sample.drop_duplicates(subset=["clone_id"])

                        tmp_a = Counter(info_sample[f"IR_VJ_1_{gene.lower()}_call"].dropna())
                        tmp_b = Counter(info_sample[f"IR_VDJ_1_{gene.lower()}_call"].dropna())

                        usage.loc[tmp_a.keys(), sample] = list(tmp_a.values())
                        usage.loc[tmp_b.keys(), sample] = list(tmp_b.values())

                    for sample in samples:
                        try:
                            chains = dict_chains[cell] if cell == "TCR" else dict_chains[cell][gene]
                            for chain in chains:
                                mask = usage.index.str.startswith(chain)
                                total = np.sum(usage.loc[mask, sample])
                                usage.loc[mask, sample] /= total
                            usage[sample] *= 100
                        except ZeroDivisionError:
                            continue

                    # Formatting output
                    usage["#Chr"] = "6"
                    usage["start"] = np.arange(1, len(usage) + 1)
                    usage["end"] = usage["start"] + 1
                    usage["ID"] = usage.index.str.replace("*", "_")

                    output_filename = f"phenotypes/usage.{cell}.{gene}.qc.functional.{weight_method}.{full}.{celltype}.{adjust}.bed.gz"
                    usage.to_csv(output_filename, sep="\t", header=True, index=False)

                    # For Plink
                    df = pd.DataFrame(index=samples, columns=["IID"])
                    df["IID"] = df.index

                    usage.set_index("ID", inplace=True)
                    df = df.join(usage.loc[:, usage.columns[4:]], how="left")
                    df.columns = df.columns.str.replace(r'(.*)/(.*)', r'\1_\2', regex=True)

                    qt = QuantileTransformer(random_state=0, output_distribution='normal')
                    df.iloc[:, 1:] = qt.fit_transform(df.iloc[:, 1:])
                    df.insert(0, "FID", "0")

                    output_filename = f"phenotypes/usage.{cell}.{gene}.qc.functional.{weight_method}.{full}.{celltype}.{adjust}.ranknorm.pheno.gz"
                    df.to_csv(output_filename, sep="\t", header=True, index=False)

