install.packages(c("ggplot2", "dplyr", "tidyr", "pheatmap", "cowplot",
                   "RColorBrewer", "ggnewscale", "ape", "factoextra", "ggfortify"),
                 repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c("limma", "EnhancedVolcano", "ComplexHeatmap",
                       "ggtree", "ggtreeExtra"))
