library(DESeq2)
library(tidyverse)
library(knitr)
library(stringr)

# Set directory and file names
base_dir <- ""
sample_list <- ""
setwd(base_dir)

args <- commandArgs(trailingOnly = TRUE)

# 10000 iterations were split into 200 parts
size_iter <- 50
split_iter <- args[1]
iter_start <- 1 + size_iter*(split_iter-1)
iter_end <- size_iter*split_iter

samples <- scan(sample_list,
                what = character())
num_sample <- length(samples)

dosage <- read.delim(dosage_file, sep = "\t")
dosage <- dosage[, c("X0", "ref", "alt", samples)]
info_genotypes <- dosage[, 1:3]
genotypes <- dosage[, 4:(3 + length(samples))]
num_genotypes <- nrow(genotypes)

is_classicalHLA <- substr(info_geno[, 1], 1, 5) %in% paste0("HLA_", c("A", "B", "C")) |
                    substr(info_geno[, 1], 1, 8) %in% paste0("HLA_", c("DRB1", "DPB1", "DPA1", "DQB1", "DQA1"))

is_classicalAA <- substr(info_geno[, 1], 1, 4) %in% paste0("AA_", c("A", "B", "C")) |
                   substr(info_geno[, 1], 1, 7) %in% paste0("AA_", c("DRB1", "DPB1", "DPA1", "DQB1", "DQA1"))

is_classical <- is_classicalHLA | is_classicalAA
index_test <- which(is_classical)

covar <- read.delim(covariate_file, sep = "\t")
rownames(covar) <- covar$sample.id
covar <- covar[samples, ]

exp <- read.delim("HLAassoc_permutation/data/sim_data.txt.gz",
                 sep = "\t")
rownames(exp) <- exp$ID
exp <- exp[paste0("iter_", seq(iter_start, iter_end)),]
exp <- exp[,samples]
num_exp <- dim(exp)[1]
covar["case"] <- covar["Status"]
include <- !is.na(covar$case) & !is.na(covar$PC1_g)
include <- include & complete.cases(t(exp))
name_dir <- file.path("HLAassoc_permutation", "results", 
                     paste(cell, split_iter, sep="."))
if (!dir.exists(name_dir)) {
  dir.create(name_dir, recursive = TRUE)
}

# Create data frame for DESeq2
i = 5
df <- cbind(covar, as.numeric(t(geno[i,])))
colnames(df)[dim(df)[2]] <- "variant"

design_formula <- "~variant+case+Age+Sex+Version+PC1_g+PC2_g"

dds <- DESeqDataSetFromMatrix(round(exp[, include]),
                              colData = df[include],
                              design = as.formula(design_formula))

dds_est <- dds
normalizationFactors(dds_est) <- matrix(1, ncol = ncol(dds_est), nrow = nrow(dds_est))
dds_est <- estimateDispersionsGeneEst(dds_est)
dispersions(dds_est) <- mcols(dds_est)$dispGeneEst

for (j in index_test) {
  name_variant <- gsub(":", "_", as.character(info_geno[j, 1]))
  name_file = file.path(name_dir, paste("permutation", cell, split_iter, name_variant, "result.txt", sep="."))
  if (file.exists(paste0(name_file))) next

  df$variant <- as.numeric(unlist(t(geno[j, ])))

  if (var(df$variant[include]) == 0) next  # Skip if no variation

  if (info_geno[j, 2] == "A") {
    dds_est$variant <- 2 - df$variant[include]
  } else if (info_geno[j, 2] == "P") {
    dds_est$variant <- df$variant[include]
  }

  # Run Wald test
  run_wald <- tryCatch({
    dds_Wald <- nbinomWaldTest(dds_est)
    results(dds_Wald, name = "variant")
  }, error = function(e) {
    message("Skipping variant due to error: ", name_variant, " - ", e$message)
    return(NULL)
  })

  write.table(run_wald, file = name_file, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)  
}
