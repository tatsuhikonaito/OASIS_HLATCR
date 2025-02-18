library(DESeq2)
library(tidyverse)
library(knitr)
library(stringr)

# Set working directory
base_dir <- ""
sample_list <- ""
dosage_file <- ""
covariate_file <- ""

setwd(base_dir)

args <- commandArgs(trailingOnly = TRUE)

cell <- args[1]
gene <- args[2]
norm <- args[3]
celltype <- args[4]
test <- args[5]
adjust <- args[6]
full <- args[7]


samples <- scan(sample_list,
                what = character())
num_sample <- length(samples)

dosage <- read.delim(dosage_file, sep = "\t")
dosage <- dosage[, c("X0", "ref", "alt", samples)]
info_genotypes <- dosage[, 1:3]
genotypes <- dosage[, 4:(3 + length(samples))]
num_genotypes <- nrow(genotypes)

is_classicalHLA <- substr(info_genotypes[, 1], 1, 5) %in% paste("HLA_", c("A", "B", "C"), sep = "") |
  substr(info_genotypes[, 1], 1, 8) %in% paste("HLA_", c("DRB1", "DPB1", "DPA1", "DQB1", "DQA1"), sep = "")
is_classicalAA <- substr(info_genotypes[, 1], 1, 4) %in% paste("AA_", c("A", "B", "C"), sep = "") |
  substr(info_genotypes[, 1], 1, 7) %in% paste("AA_", c("DRB1", "DPB1", "DPA1", "DQB1", "DQA1"), sep = "")
is_classical <- is_classicalHLA | is_classicalAA
index_test <- which(is_classical)

covar <- read.delim(covariate_file, sep = "\t")
rownames(covar) <- as.vector(unlist(covar[,"sample.id"]))
covar <- covar[samples,]

prefix <- paste(cell, gene, "qc.functional", norm, full, celltype, adjust, sep = ".")
exp <- read.delim(file.path("phenotypes", paste("usage", prefix, "bed.gz", sep = ".")), sep = "\t")
rownames(exp) <- exp$ID
info_tcr <- exp[, 4]
exp <- exp[, samples]
num_exp <- nrow(exp)
include <- !is.na(covar$Status) & !is.na(covar$PC1_g)
include <- include & complete.cases(t(exp))
covar["case"] <- covar["Status"]
num_case <- sum(covar$case == 1)

include <- !is.na(covar$case) & !is.na(covar$PC1_g)
name_dir <- file.path(paste("results_HLAassoc", sep = "."), paste("usage", prefix, test, sep = "."))
if (!dir.exists(name_dir)) {
  dir.create(name_dir)
}

df <- cbind(covar, c(rep(0, num_case/2), rep(1, num_case/2), rep(0, (num_sample - num_case)/2), rep(1, (num_sample - num_case)/2)))
colnames(df)[ncol(df)] <- "variant"

if (test == "qtl") {
  design_formula <- "~variant+case+Age+Sex+Version+PC1_g+PC2_g"
} else if (test == "int") {
  design_formula <- "~case+variant+case:variant+Age+Sex+Version+PC1_g+PC2_g"
}

dds <- DESeqDataSetFromMatrix(round(exp[, include]),
                              colData = df[include,],
                              design = as.formula(design_formula))

dds_est <- dds
normalizationFactors(dds_est) <- matrix(rep(1, ncol(dds_est) * nrow(dds_est)),
                                         ncol = ncol(dds_est), nrow = nrow(dds_est))
dds_est <- estimateDispersionsGeneEst(dds_est)
dispersions(dds_est) <- mcols(dds_est)$dispGeneEst

for (j in index_test) {
  name_variant <- gsub(":", "_", as.character(info_genotypes[j, 1]))
  name_file <- paste("usage.deseq2", prefix, test, name_variant, "result.txt", sep = ".")
  if (file.exists(paste0(name_file, ".gz"))) next
  
  df <- cbind(covar, t(genotypes[j,]))
  colnames(df)[ncol(df)] <- "variant"
  df$variant <- as.numeric(unlist(t(genotypes[j, ])))

  if (var(df$variant[include]) == 0) next  # Skip if no variation

  if (info_genotypes[j, 2] == "A") {
    dds_est$variant <- 2 - df$variant[include]
  } else if (info_genotypes[j, 2] == "P") {
    dds_est$variant <- df$variant[include]
  }
  run_wald <- try(nbinomWaldTest(dds_est), silent = FALSE)
  if (class(run_wald) == "try-error") {
    next
  }
  if (test == "qtl") {
    res <- results(dds_est, name = "variant")
  } else if (test == "int") {
    res <- results(dds_est, name = "case.variant")
  }
  write.table(res, file = file.path(name_dir, name_file),
              row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
  system(paste("gzip -f", file.path(name_dir, name_file)))
}
