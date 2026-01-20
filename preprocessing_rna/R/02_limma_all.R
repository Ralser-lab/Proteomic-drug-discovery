#########################
# 0. Packages & Paths
#########################

# install.packages("rstudioapi")
# install.packages("data.table")
# install.packages("argparse")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

# BiocManager::install("edgeR")
# BiocManager::install("limma")

library(rstudioapi)
library(data.table)
library(argparse)
library(edgeR)
library(limma)

# Parse command line arguments
parser <- ArgumentParser(description = 'Limma differential expression analysis')
parser$add_argument(
  '--norm',
  type = 'character',
  default = 'none',
  help = 'Normalization method: "none", "tmm", or "voom"'
)
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
# 1. Extract data
#########################

script_dir  <- dirname(get_script_path())
data_dir    <- args$data_dir

data_df <- as.data.frame(fread(
  file.path(script_dir, '..', data_dir, "raw_counts_filtered_TMM_logCPM.csv")
))

# Format data_df index
row.names(data_df) <- data_df$V1
data_df$V1 <- NULL 
data_df <- t(data_df)

meta_df <- as.data.frame(fread(
  file.path(script_dir, '..', data_dir, "metadata_limma_design.csv"),
  header = TRUE
))

#Format meta_df index
row.names(meta_df) <- meta_df$V1
meta_df$V1 <- NULL

contrast_df <- as.data.frame(fread(
  file.path(script_dir, '..', data_dir, "contrasts.csv"),
  header = TRUE
))


#########################
# 1. Create design and contrast matrix
#########################

condition_ <- factor(meta_df$condition)

design_matrix <- model.matrix(~0 + condition_)
colnames(design_matrix) <- levels(condition_)

contrast_matrix <- makeContrasts(
  contrasts = contrast_df$contrast,
  levels    = design_matrix
)

#########################
# 5. Limma: branch by --norm
#########################

cat('Performing dummy regression...')
fit1 <- lmFit(data_df, design_matrix)

cat('Fitting contrasts...')
fit2 <- contrasts.fit(fit1, contrast_matrix)

cat('Estimating variance (sample var borrowed from genes var)...')
fit2 <- eBayes(fit2, trend = TRUE)

print('Writing limma obj to disk (.rds)...')
saveRDS(fit2, 
        file.path(script_dir, "..", "data", "tmm_limma_fit2_ebayes.rds")
)

