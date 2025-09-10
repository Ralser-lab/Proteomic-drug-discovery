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

#Set absolute path to /data folder in clone github repository
abspath = '/Users/shaon/Desktop/PROTACS/github_deposition/data'
setwd(abspath)

plot_df = read.csv2('analog_dataout.csv', sep = ',', dec = '.')

metadata = plot_df[ncol(plot_df)]

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

setwd(paste(abspath, '/../figures/', sep = ''))

pdf('heatmap_analogues_weights_predictions.pdf', width = 4, height = 3)

pheatmap(t(de_matrix), color = color, cluster_rows = TRUE, cluster_cols = FALSE,
         annotation_col = annotation_col, 
         annotation_colors = annotation_colors, border_color = 'black')

dev.off()


