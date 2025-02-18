#!/bin/bash
# Possible variable values:
#   pheno: "usage", "pc", "cdr3_len", "mu_freq"
#   cell: "BCR" or "TCR"
#   If cell == "BCR", then celltype can be one of:
#         "B", "B_Activated", "B_Intermediate", "B_Memory", "B_Naive", "B_Naive1", "B_Naive2", "PB"
#   If cell == "TCR", then celltype can be one of:
#         "CD4T", "CD8T", "Other_T", "CD4_CTL", "CD4_Memory", "CD4_Naive", "CD8_CTL", "CD8_Memory", "CD8_Naive", "MAIT", "Pro_T", "Treg", "CD4_TCM", "CD4_TEM", "CD8_TCM", "CD8_TEM"
#   gene: "V", "D", "J"
#   If cell == "BCR", then chain can be one of:
#         "H", "K", "L"
#   If cell == "TCR", then chain can be one of:
#         "A", "B"

MEM=""
THREADS=""

PLINK2="plink2 --memory ${MEM}000 --threads ${THREADS} "

GENO=""
COVAR=""

pheno="$1"
cell="$2"
shift 2  # Move past $1 (pheno) and $2 (cell)

# Fixed parameters
weight_method="raw_cdr3"
adjust="norm"

case "$pheno" in
    "usage")
        celltype="$1"
        gene="$2"
        full="full"
        PHENO="phenotypes/usage.${cell}.${gene}.qc.functional.${weight_method}.${full}.${celltype}.${adjust}.ranknorm.pheno.gz"
        OUT="results_GWAS/usage.${cell}.${gene}.qc.functional.${weight_method}.${full}.${celltype}.${adjust}"
        ;;
    "pc")
        celltype="$1"
        gene="$2"
        chain="$3"
        full="full"
        PHENO="phenotypes/pc.${cell}.${gene}.qc.functional.${weight_method}.${full}.${celltype}.${adjust}.${chain}.ranknorm.pheno.gz"
        OUT="results_GWAS/pc.${cell}.${gene}.qc.functional.${weight_method}.${full}.${celltype}.${adjust}.${chain}"
        ;;
    "cdr3_len"|"mu_freq")
        PHENO="phenotypes/${pheno}.${cell}.ranknorm.pheno.gz"
        OUT="results_GWAS/${pheno}.${cell}"
        ;;
    *)
        echo "Error: Invalid phenotype type '${pheno}'."
        exit 1
        ;;
esac


${PLINK2} --pfile "${GENO}" \
    --linear hide-covar \
    --pheno "${PHENO}" \
    --covar "${COVAR}" \
    --covar-variance-standardize \
    --out "${OUT}"
