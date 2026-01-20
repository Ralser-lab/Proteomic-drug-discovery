#########################
# 0. Packages & Paths
#########################

library(rstudioapi)
library(data.table)
library(argparse)
library(ggplot2)
library(limma)

get_script_path <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    return(rstudioapi::getActiveDocumentContext()$path)
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args)
  if (length(file_arg) > 0)
    return(normalizePath(sub("--file=", "", args[file_arg])))
  if (!is.null(sys.frames()[[1]]$ofile))
    return(normalizePath(sys.frames()[[1]]$ofile))
  return(NULL)
}

parser <- ArgumentParser(description = 'Generate volcano plots (one per contrast)')
parser$add_argument(
  "--data-dir", type = "character", default = "data",
  help = 'Directory name for reading data files (default: "data")'
)
args <- parser$parse_args()

#########################
# 1. Extract data
#########################

script_dir <- dirname(get_script_path())
data_dir   <- args$data_dir

fit2 <- readRDS(file.path(script_dir, "..", data_dir, "tmm_limma_fit2_ebayes.rds"))

# where to save plots
out_dir <- file.path(script_dir, "..", "plots", "volcano_by_contrast")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#########################
# 2. Loop over contrasts
#########################

contrast_names <- colnames(fit2$coefficients)

if (is.null(contrast_names) || length(contrast_names) == 0) {
  stop("No contrasts found in fit2$coefficients. Did you run contrasts.fit()?")
}

for (co in contrast_names) {
  message("Processing contrast: ", co)
  
  tt <- topTable(fit2, coef = co, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  
  # volcano data with direction
  volcano_df <- data.frame(
    gene = tt$gene,
    logFC = tt$logFC,
    neglog10FDR = -log10(tt$adj.P.Val),
    direction = ifelse(
      tt$adj.P.Val < 0.05 & tt$logFC >= log2(1.5), "Up",
      ifelse(tt$adj.P.Val < 0.05 & tt$logFC <= -log2(1.5), "Down", "NS")
    )
  )
  
  p <- ggplot(volcano_df, aes(x = logFC, y = neglog10FDR, color = direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(
      values = c(
        "Up" = "darkred",
        "Down" = "darkblue",
        "NS" = "grey70"
      )
    ) +
    labs(
      x = "log2 fold change",
      y = "-log10(adj. p-value)",
      title = co,
      color = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    ) +
    coord_cartesian(xlim = c(-6, 6), ylim = c(0, 15))
  
  # safe filename per contrast
  safe_co <- gsub("[^A-Za-z0-9]+", "_", co)
  file_png <- file.path(out_dir, paste0("volcano_", safe_co, ".png"))
  
  ggsave(
    filename = file_png,
    plot = p,
    width = 4.5,
    height = 4.0,
    dpi = 300
  )
  
  message("  Saved: ", file_png)
}

