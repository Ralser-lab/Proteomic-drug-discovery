#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#'
#' Script Name: protacs_06_limma_drug_10uM.R
#' Description: Differential expression analysis of PROTAC proteomes at 10 µM
#'              using limma (drug vs DMSO contrasts), with export of matrices
#'
#' Author: Shaon Basu
#' Date: 2025-09-16
#'
#' Inputs
#' ------
#' - data/SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv
#' - data/SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_metaforlimma_240611a.tsv
#'
#' Outputs
#' -------
#' - data/Drug_LFCxPval_250305a.csv         
#' - data/Drug_LFCxadjPval_250305a.csv      
#'
#' Requirements
#' ------------
#' R >= 4.2
#' Packages: BiocManager, limma, EnhancedVolcano, ggplot2, dplyr, tidyr,
#'           pheatmap, grid, cowplot
#'

# Define required CRAN and Bioconductor packages
cran_packages <- c("ggplot2", "dplyr", "tidyr", "pheatmap", "grid", "cowplot")
bioc_packages <- c("limma", "EnhancedVolcano")

# Function to install missing CRAN packages
install_if_missing <- function(packages) {
  installed <- rownames(installed.packages()) # Get installed packages
  missing <- packages[!(packages %in% installed)] # Check missing ones
  if (length(missing) > 0) {
    message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org/")
  }
}

# Ensure Bioconductor is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Install missing CRAN packages
install_if_missing(cran_packages)

# Install missing Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Load all required libraries
required_packages <- c(cran_packages, bioc_packages)
lapply(required_packages, library, character.only = TRUE)

# Set up relative paths
data_dir <- file.path(getwd(), "data")
fig_dir  <- file.path(getwd(), "figures")
setwd(data_dir)

# Load metadata and expression data (limma formatted proteomes) from HBD screen
metadata <- read.csv2('SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_metaforlimma_240611a.tsv', sep = ',', dec = '.', row.names = 1)

expressionmatrix <- read.csv2('SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv', 
                              sep = ',', dec = '.', row.names = 1, header = 1,
                              check.names = FALSE)

# Format expression and metadata for perfomring limma
colnames(metadata) <- c('Compound_', 'Concentration_', 'Plate_', 'Check_', 'SpRep_', 'SpBatch_', 'SpPosition_', 'Annotation_', "IC50_Glu_", "IC50_Gal_",
                        'E3_ligase_', 'Target_', 'Drug_Type_', 'Binned_Ligase_', 'Binned_Target_', 'Cluster_', 'Dend_','Cluster2_')
metadata[is.na(metadata)] <- 0
metadata[] <- lapply(metadata, function(col) {
  if (is.factor(col)) {
    as.factor(gsub("-", "_", as.character(col)))
  } else if (is.character(col)) {
    gsub("-", "_", col)
  } else {
    col
  }
})

# Subset 10 micromolar data
metadata_10_DMSO <- subset(metadata, metadata$Concentration_ == 0.1 | metadata$Concentration_ == 1 |
                         metadata$Concentration_ == 10 | metadata$Concentration_ == 0)
metadata_10_DMSO <- subset(metadata, metadata$Concentration_ == 10 | metadata$Concentration_ == 0)
expression_10_DMSO <- expressionmatrix[, colnames(expressionmatrix) %in% rownames(metadata_10_DMSO)]
z <- t(table(metadata_10_DMSO$Compound_, metadata_10_DMSO$Cluster2_))

# Format metadata for drug-wise contrasts
metadata_10_DMSO$Cluster2_ <- as.factor(metadata_10_DMSO$Cluster2_)
metadata_10_DMSO$Plate_ <- as.factor(metadata_10_DMSO$Plate_)
indices <- metadata_10_DMSO$Cluster2_ == 'DMSO'  # Find indices where Cluster2_ is 'DMSO'
metadata_10_DMSO$Annotation_[indices] <- 'DMSO' 
metadata_10_DMSO$Compound_[indices] <- 'DMSO' 
metadata_10_DMSO$Annotation_ <- as.factor(metadata_10_DMSO$Annotation_)
metadata_10_DMSO$Compound_ <- as.factor(metadata_10_DMSO$Compound_)

# Format annotation to be compatible with R dataframe
levels(metadata_10_DMSO$Annotation_) <- make.names(levels(metadata_10_DMSO$Annotation_))
design <- model.matrix(~0+Compound_,metadata_10_DMSO)

# Fit limma model for fold change modelling
fit <- lmFit(expression_10_DMSO, design)

# Build drug contrasts vs DMSO (each drug in HBD library)
drugs <- levels(metadata_10_DMSO$Compound)[levels(metadata_10_DMSO$Compound) != "DMSO"]
drugs <- paste0("Compound_", drugs)
contrast_formulas <- setNames(sapply(drugs, function(drug) paste0(drug, " - Compound_DMSO"), USE.NAMES = FALSE), 
                              sapply(drugs, function(drug) paste0(drug, "vsCompound_DMSO")))
contrasts <- makeContrasts(contrasts = as.list(contrast_formulas), levels = design)

# Perform limma test to calculate probabilities (significance)
fit2 <- contrasts.fit(fit, contrasts)
fit3 <- eBayes(fit2, trend = TRUE)

# Extract summary statistics and probabilities
pval <- fit3$p.value
adjpval <- apply(pval, 2, function(x) p.adjust(x, method = "BH"))
logFC_matrix2 <- fit3$coefficients
limma_matrix2 <- -log10(pval) * sign(logFC_matrix2)
limma_matrix3 <- -log10(adjpval) * sign(logFC_matrix2)

# Export differential expression data from limma (as matrices)
write.csv2(limma_matrix2, 'Drug_LFCxPval_250305a.csv')
write.csv2(limma_matrix3, 'Drug_LFCxadjPval_250305a.csv')


