#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#'
#' Script Name: protacs_19_pheatmap_topweights_all.R
#' Description: 
#'   Generate heatmaps of final model outputs (probabilities, predictions, 
#'   and feature weights) from proteome-based HBD screens. This script 
#'   processes model-derived tables, builds annotations, and visualizes 
#'   feature expression and prediction profiles as clustered heatmaps.
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
#' - figures/protacs_19_pheatmap_topweights_all.pdf
#'   (Clustered heatmap of feature matrix with annotation tracks)
#'
#' Requirements
#' ------------
#' R >= 4.2
#' CRAN packages: ggplot2, RColorBrewer, pheatmap
#' Bioconductor packages: ComplexHeatmap
#'
#'

# Load required CRAN and Bioconductor packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
})

# Set up relative paths
data_dir <- file.path(getwd(), "data")
fig_dir  <- file.path(getwd(), "figures")
setwd(data_dir)

# Read in xgb finalmodel input/output related data
plot_df = read.csv2('predictions_shap.csv', sep = ',', dec = '.')
plot_df$Actual <- as.factor(plot_df$Actual)
plot_df$Predicted <- as.factor(plot_df$Predicted)
plot_df$Drug <- as.factor(plot_df$Drug)
plot_df$Cluster <- as.factor(plot_df$Cluster)

# Extract metadata for R heatmap
metadata = plot_df[1:6]

# Create annotations for pheatmap
annotation_col <- data.frame(
                           Probability = metadata$Probability,
                           Toxicity = metadata$Actual,
                             Drug = metadata$Drug)

row.names(annotation_col) <- plot_df$X
de_matrix = plot_df[7:ncol(plot_df)]
row.names(de_matrix) = metadata$X
de_matrix = as.matrix(de_matrix)
annotation_row <- data.frame(ID = c('ComplexI'))
binary_color = c('blue', 'red')
proba_color = c('black', 'green')
predicted_levels <- levels(metadata$Predicted)
actual_levels <- levels(metadata$Actual)
predicted_colors <- setNames(c(binary_color), predicted_levels)
actual_colors <- setNames(c(binary_color), actual_levels)
drug_levels <- levels(metadata$Drug)
spectral_colors <- colorRampPalette(brewer.pal(11, "Spectral"))(length(drug_levels))
drug_colors <- setNames(spectral_colors, drug_levels)
tox_colors <- color <- colorRampPalette((c(proba_color)))(100)

annotation_colors <- list(
  Toxicity = actual_colors,
  Prediction = predicted_colors,
  Probability = tox_colors,
  Drug = drug_colors
)

color <- colorRampPalette((c('darkblue','darkblue','darkblue','blue','blue','blue','blue','blue',"white", "red",'red','darkred')))(213)
ordered <- de_matrix[order(de_matrix[,1],decreasing=TRUE),]

# Plot and save R pheatmap
setwd(fig_dir)
pdf('protacs_19_pheatmap_topweights_all.pdf', width = 10, height = 5)
pheatmap(t(de_matrix), color = color, scale = 'none', cluster_rows = TRUE, cluster_cols = TRUE, 
        annotation_col = annotation_col, 
        show_colnames = FALSE,
         annotation_colors = annotation_colors, border_color = 'black')
dev.off()


