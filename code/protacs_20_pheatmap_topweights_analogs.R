#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#'
#' Script Name: protacs_20_pheatmap_topweights_analogs.R
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
#' - figures/protacs_20_pheatmap_topweights_analogs.pdf
#'   (Clustered heatmap of feature matrix with annotation tracks)
#'
#' Requirements
#' ------------
#' R >= 4.2
#' Packages: ggplot2, RColorBrewer, pheatmap
#' Bioconductor packages: ComplexHeatmap
#'
#'

# Load required CRAN and Bioconductor packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
  library(ComplexHeatmap)
})

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
de_matrix = plot_df[c('NDUFA5','CYC1','NDUFA4')]
row.names(de_matrix) = plot_df$analogue
de_matrix = as.matrix(de_matrix)
proba_color = c('black', 'green')
tox_colors <- color <- colorRampPalette((c(proba_color)))(100)
annotation_colors <- list(
  Probability = tox_colors
)
# Asymmetric ramp: blue at -0.5, white at 0, red at 1
# Proportional split: 33 colours for [-0.5, 0], 67 colours for [0, 1]
n_neg   <- 33
n_pos   <- 67
color   <- c(colorRampPalette(c("blue",  "white"))(n_neg),
             colorRampPalette(c("white", "red"  ))(n_pos))
breaks  <- c(seq(-0.5, 0, length.out = n_neg + 1),
             seq( 0,   1, length.out = n_pos + 1)[-1])
ordered <- de_matrix[order(de_matrix[,1],decreasing=TRUE),]

# Plot and save R pheatmap on series 15, compound 1analogues subset
setwd(fig_dir)
pdf('protacs_20_pheatmap_topweights_analogs.pdf', width = 4, height = 3)
pheatmap(t(de_matrix), color = color, breaks = breaks, cluster_rows = TRUE, cluster_cols = FALSE,
         annotation_col = annotation_col, 
         annotation_colors = annotation_colors, border_color = 'black')
dev.off()


