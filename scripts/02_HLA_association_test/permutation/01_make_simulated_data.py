#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm

# Set base directory
base_dir = ""
os.chdir(base_dir)

seed_value = ""
np.random.seed(seed_value)

# Parameters for a negative binomial regression were determined as median values of those fitted for actual TCR counts in representative settings.
adjust = "norm"
cell = "TCR"
gene = "V"
weight_method = "weighted"
full = "full"
celltype = "all"

# File path construction
input_file = f"phenotypes/{cell}.{gene}.qc.functional.{weight_method}.{full}.{celltype}.{adjust}.bed.gz"

# Load data
usage = pd.read_csv(input_file, index_col=3, sep="\t").iloc[:, 3:]

# Count zero values in each row
num_sample = usage.shape[1]
usage["num_0"] = (usage == 0).sum(axis=1)

# Filter genes with at least 80% nonzero values
pass_usage = usage["num_0"] < num_sample * 0.2

df = pd.DataFrame(index = exp.index[pass_usage],
                  columns = ["n", "p"])

# Fit Negative Binomial model
for i in df.index:
    try:
        tmp = usage.loc[i].iloc[3:]
        x = np.ones_like(tmp)
        res = sm.NegativeBinomial(tmp, x).fit(start_params=[1, 1], disp=False)
        
        mu = np.usage(res.params[0])
        p = 1 / (1 + np.usage(res.params[0]) * res.params[1])
        n = np.usage(res.params[0]) * p / (1 - p)
        
        df.loc[i] = [n, p]

    except Exception as e:
        print(f"Warning: Model fitting failed for index {i}: {e}")

n_median, p_median = df.median()

n_iter = 10000

# Generate simulated data
usage_sim = pd.DataFrame(
    np.random.negative_binomial(n_median, p_median, (n_iter, num_sample)),
    index=[f"iter_{i}" for i in range(1, n_iter + 1)],
    columns=usage.columns[:num_sample]
)
usage_sim.index.name = "ID"

# Output file paths
output_dir = "HLAassoc_permutation/data"
os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists
usage_sim.to_csv(f"{output_dir}/sim_data.txt.gz", sep="\t", header=True, index=True)
