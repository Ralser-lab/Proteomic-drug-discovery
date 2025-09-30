#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#'
#' Script Name: protacs_19_pheatmap_topweights_analogs.R
#' Description: 
#'   Visualize proteomic feature signatures and final model-derived scores 
#'   for HBD compound analogues (Series 15, Compound 1 chemical series). 
#'   This script loads analogue-level feature means, formats metadata into 
#'   annotations, and generates clustered heatmaps to compare analogue 
#'   signatures and toxicity scores.
#'
#' Author: Shaon Basu
#' Date: 2025-09-30
#'
#' Inputs
#' ------
#' - data/Rplot_Figure4.csv
#'   (Model predictions, probabilities, and SHAP-derived feature matrix)
#'
#' Outputs
#' -------
#' - figures/protacs_19_pheatmap_topweights_analogs.pdf
#'   (Clustered heatmap of feature matrix with annotation tracks)
#'
#' Requirements
#' ------------
#' R >= 4.2
#' Packages: ggplot2, RColorBrewer, pheatmap
#' Bioconductor packages: ComplexHeatmap
#'
#'

# Define required CRAN and Bioconductor packages
cran_packages <- c("ggplot2", "RColorBrewer", "pheatmap")
bioc_packages <- c("ComplexHeatmap")

# Function to install missing CRAN packages
install_if_missing <- function(packages) {
  installed <- rownames(installed.packages()) # Get installed packages
  missing <- packages[!(packages %in% installed)] # Check for missing ones
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

message("All required packages are installed and loaded successfully!")

# Set up relative paths
data_dir <- file.path(getwd(), "data")
fig_dir  <- file.path(getwd(), "figures")
setwd(data_dir)

# Read in analogue related input/ouptus and related data
plot_df = read.csv2('analog_dataout.csv', sep = ',', dec = '.')

# Extract metadata
metadata = plot_df[ncol(plot_df)]

# Format annotation data for R pheatmap
annotation_col <- data.frame(
                           Probability = metadata$probability)
row.names(annotation_col) <- plot_df$analogue
de_matrix = plot_df[c('NDUFA5','NDUFA4', 'NDUFB10', 'NDUFA13')]
row.names(de_matrix) = plot_df$analogue
de_matrix = as.matrix(de_matrix)
proba_color = c('blue', 'red')
tox_colors <- color <- colorRampPalette((c(proba_color)))(100)
annotation_colors <- list(
  Probability = tox_colors
)
color <- colorRampPalette((c('blue',"white", "red")))(100)
ordered <- de_matrix[order(de_matrix[,1],decreasing=TRUE),]

# Plot and save R pheatmap on series 15, compound 1analogues subset
setwd(fig_dir)
pdf('protacs_19_pheatmap_topweights_analogs.pdf', width = 4, height = 3)
pheatmap(t(de_matrix), color = color, cluster_rows = TRUE, cluster_cols = FALSE,
         annotation_col = annotation_col, 
         annotation_colors = annotation_colors, border_color = 'black')
dev.off()


