#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
import numpy as np

# Set directory and file names
base_dir = ""
HLA_variant_list_file = ""
output_file = "p_threshold.txt"
os.chdir(base_dir)

n_iter = 10000

variants = np.loadtxt(HLA_variant_list_file)

results = pd.read_csv("permutation/results/permutation_summary.TCR.1_200.result.txt.gz",
                     sep="\t", header=0)

results["variant"] = results.FN.str.split(".", expand=True)[4]
results = results.loc[results.variant.isin(variants)]
variants = sorted(set(results.variant))

df = pd.DataFrame(index=[f"iter_{i}" for i in range(1, n_iter+1)], columns=variants)

for variant in df.columns:
    df[variant] = np.array(results.loc[results.variant == variant, "P"]).flatten()
    
pmins = np.array(np.min(df, axis=1))
p_thr = np.percentile(pmins, 5)

with open(output_file, "w") as f:
    f.write(str(p_thr) + "\n")

print(f"p_thr saved to {output_file}: {p_thr}")
