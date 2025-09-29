# Define required CRAN and Bioconductor packages
cran_packages <- c("ggplot2", "RColorBrewer", "pheatmap")
bioc_packages <- c("ComplexHeatmap")

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

message("All required packages are installed and loaded successfully!")

#Set absolute path to /data folder in clone github repository
abspath = '/Users/shaon/Desktop/PROTACS/github_deposition/data'
setwd(abspath)

plot_df = read.csv2('Rplot_Figure4.csv', sep = ',', dec = '.')

plot_df$Actual <- as.factor(plot_df$Actual)
plot_df$Predicted <- as.factor(plot_df$Predicted)
plot_df$Drug <- as.factor(plot_df$Drug)
plot_df$Cluster <- as.factor(plot_df$Cluster)

metadata = plot_df[1:6]

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
proba_color = c('blue', 'red')

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

setwd(paste(abspath, '/../figures/', sep = ''))

pdf('heatmap_all_weights_predictions.pdf', width = 10, height = 4)

pheatmap(t(de_matrix), color = color, scale = 'none', cluster_rows = TRUE, cluster_cols = TRUE, 
        # annotation_col = annotation_col, 
        # labels_col = NA,
         annotation_colors = annotation_colors, border_color = 'black')

dev.off()


