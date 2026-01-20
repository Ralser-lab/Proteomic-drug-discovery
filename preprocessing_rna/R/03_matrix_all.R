#########################
# 0. Packages & Paths
#########################

# install.packages("rstudioapi")
# install.packages("data.table")
# install.packages("argparse")

library(rstudioapi)
library(data.table)
library(argparse)

# Parse command line arguments
parser <- ArgumentParser(description = 'Extract matrices from limma fit object')
parser$add_argument('--norm', type = 'character', default = 'tmm',
                   help = 'Normalization method: "none" or "tmm" or "voom')
parser$add_argument(
  "--data-dir", type = "character", default = "data",
  help = 'Directory name for reading/writing data files (default: "data")'
)
args <- parser$parse_args()

get_script_path <- function() {
  # RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    return(rstudioapi::getActiveDocumentContext()$path)
  }
  
  # Rscript
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args)
  if (length(file_arg) > 0)
    return(normalizePath(sub("--file=", "", args[file_arg])))
  
  # source()
  if (!is.null(sys.frames()[[1]]$ofile))
    return(normalizePath(sys.frames()[[1]]$ofile))
  
  # fallback
  return(NULL)
}

#########################
# 1. Extract and transform dataframes from limma fit object
#########################

script_dir  <- dirname(get_script_path())
data_dir    <- args$data_dir

# Determine suffix based on normalization method
if (args$norm == "tmm") {
  norm_suffix <- "tmm_"
  norm_label  <- "TMM_"
} else if (args$norm == "voom") {
  norm_suffix <- "voom_"
  norm_label  <- "VOOM_"
} else {
  norm_suffix <- ""
  norm_label  <- ""
}

fit2 <- readRDS(file.path(script_dir
                         , '..', data_dir, 
                         paste0(norm_suffix,'tmm_limma_fit2_ebayes.rds')))

lfc_matrix <- fit2$coefficients

pval_matrix <- fit2$p.value 

bh_matrix <- apply(pval_matrix, 2, p.adjust, method = "BH")

#########################
# 2. Load matrices as .csvs and write
#########################

print('Write lfc matrix...')
write.csv(
  t(lfc_matrix),
  file.path(script_dir, "..", data_dir, 
           paste0("raw_counts_filtered_", norm_label, "logCPM_limma_dmso_lfc.csv"))
)

print('Write pval matrix...')
write.csv(
  t(pval_matrix),
  file.path(script_dir, "..", data_dir, 
           paste0("raw_counts_filtered_", norm_label, "logCPM_limma_dmso_pval.csv"))
)

print('Write bh matrix...')
write.csv(
  t(bh_matrix),
  file.path(script_dir, "..", data_dir, 
           paste0("raw_counts_filtered_", norm_label, "logCPM_limma_dmso_bh.csv"))
)

