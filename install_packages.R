if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")

cran_packages <- list(
    ggplot2    = "3.5.2",
    dplyr      = "1.1.4",
    tidyr      = "1.3.1",
    pheatmap   = "1.0.13",
    cowplot    = "1.2.0",
    RColorBrewer = "1.1-3",
    ggnewscale = "0.5.2",
    ape        = "5.8-1",
    factoextra = "1.0.7",
    ggfortify  = "0.4.18"
)

for (pkg in names(cran_packages)) {
    remotes::install_version(pkg, version = cran_packages[[pkg]],
                             repos = "https://cloud.r-project.org",
                             upgrade = "never")
}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(version = "3.18", ask = FALSE, update = FALSE)

BiocManager::install(c("limma",            # 3.58.1
                       "EnhancedVolcano",  # 1.20.0
                       "ComplexHeatmap",   # 2.18.0
                       "ggtree",           # 3.10.1
                       "ggtreeExtra"),     # 1.12.0
                     ask = FALSE, update = FALSE)
