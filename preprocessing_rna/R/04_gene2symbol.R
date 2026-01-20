#########################
# 0. Packages & Paths
#########################

# install.packages("rstudioapi")
# install.packages("data.table")
# install.packages("argparse")
# install.packages("ggplot2")

library(rstudioapi)
library(data.table)
library(argparse)

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

# Parse command line arguments
parser <- ArgumentParser(description = 'Map limma models to gene symbols')
parser$add_argument(
  "--norm", type = "character", default = "tmm",
  help = 'Which normalization to process: "none", "tmm"'
)
parser$add_argument(
  "--data-dir", type = "character", default = "data",
  help = 'Directory name for reading/writing data files (default: "data")'
)
args <- parser$parse_args()

#########################
# 1. Script dir
#########################

script_dir  <- dirname(get_script_path())

###########################################################
# 1. Extract fit2 object & mart file  *TMM*
###########################################################
if (args$norm == "tmm" || args$norm == "all") {

  script_dir  <- dirname(get_script_path())
  data_dir    <- args$data_dir

  fit2 <- readRDS(file.path(script_dir
                           , '..', data_dir, 'tmm_limma_fit2_ebayes.rds'))
  
  mart <- read.delim(file.path(script_dir
                               , '..', data_dir, 'mart_GRCh38.p14.txt'), 
                     stringsAsFactors = FALSE)
  head(mart)
  colnames(mart)
  
  #########################
  # 2. Remap genes -> symbol  *TMM*
  #########################
  
  # Get limma object gene ids
  # (use the coefficients matrix, not fit2 itself)
  gene_ids <- rownames(fit2$coefficients)
  
  # Match order of fit2 genes to mart rows
  idx <- match(gene_ids, mart$Gene.stable.ID) 
  
  # Get corresponding gene symbols
  gene_symbols <- mart$Gene.name[idx]  
  
  # Replace NAs with the original gene ID so you don't lose them
  gene_symbols[is.na(gene_symbols)] <- gene_ids[is.na(gene_symbols)]
  
  # For duplicated symbols, keep the original ID instead 
  dup <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)
  gene_symbols[dup] <- gene_ids[dup]
  # now gene_symbols is mostly symbols, but duplicates have been reverted to IDs
  
  # Set rownames for all the main matrices in fit2
  slots_to_rename <- c("coefficients", "p.value", "t", "stdev.unscaled", "lods")
  
  for (s in slots_to_rename) {
    if (!is.null(fit2[[s]])) {
      rownames(fit2[[s]]) <- gene_symbols
    }
  }
  
  if (!is.null(fit2$Amean)) {
    names(fit2$Amean) <- gene_symbols
  }
  
  # Also store them in fit2$genes (keep IDs as rownames, symbols as a column)
  fit2$genes <- data.frame(
    ID     = gene_ids,
    Symbol = gene_symbols,
    stringsAsFactors = FALSE,
    row.names = gene_ids  # original IDs are unique, so no duplicate row.names error
  )
  
  #########################
  # 3. Extract dataframes  *TMM*
  #########################
  
  lfc_matrix <- fit2$coefficients
  
  pval_matrix <- fit2$p.value 
  
  bh_matrix <- apply(pval_matrix, 2, p.adjust, method = "BH")
  rownames(bh_matrix) <- rownames(pval_matrix)
  
  lfc_bh_matrix <- -log10(bh_matrix) * sign(lfc_matrix)
  rownames(lfc_bh_matrix) <- rownames(lfc_matrix)
  
  #########################
  # 4. Transform and load matrices as .csvs and write  *TMM*
  #########################
  
  print('Write lfc matrix...')
  write.csv(
    t(lfc_matrix),
    file.path(script_dir, "..", data_dir,
              "raw_counts_filtered_TMM_logCPM_limma_dmso_lfc_symbols.csv")
  )
  
  print('Write pval matrix...')
  write.csv(
    t(pval_matrix),
    file.path(script_dir, "..", data_dir,
              "raw_counts_filtered_TMM_logCPM_limma_dmso_pval_symbols.csv")
  )
  
  print('Write bh matrix...')
  write.csv(
    t(bh_matrix),
    file.path(script_dir, "..", data_dir,
              "raw_counts_filtered_TMM_logCPM_limma_dmso_bh_symbols.csv")
  )

} # end if TMM / all


###########################################################
# 1b. Extract fit2 object & mart file  *no TMM*
###########################################################
if (args$norm == "none" || args$norm == "all") {

  script_dir  <- dirname(get_script_path())
  data_dir    <- args$data_dir
  
  fit2 <- readRDS(file.path(script_dir
                           , '..', data_dir, 'limma_fit2_ebayes.rds'))
  
  mart <- read.delim(file.path(script_dir
                               , '..', data_dir, 'mart_GRCh38.p14.txt'), 
                     stringsAsFactors = FALSE)
  head(mart)
  colnames(mart)
  
  #########################
  # 2b. Remap genes -> symbol *no TMM*
  #########################
  
  # Get limma object gene ids
  # (use the coefficients matrix, not fit2 itself)
  gene_ids <- rownames(fit2$coefficients)
  
  # Match order of fit2 genes to mart rows
  idx <- match(gene_ids, mart$Gene.stable.ID) 
  
  # Get corresponding gene symbols
  gene_symbols <- mart$Gene.name[idx]  
  
  # Replace NAs with the original gene ID so you don't lose them
  gene_symbols[is.na(gene_symbols)] <- gene_ids[is.na(gene_symbols)]
  
  # For duplicated symbols, keep the original ID instead 
  dup <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)
  gene_symbols[dup] <- gene_ids[dup]
  # now gene_symbols is mostly symbols, but duplicates have been reverted to IDs
  
  # Set rownames for all the main matrices in fit2
  slots_to_rename <- c("coefficients", "p.value", "t", "stdev.unscaled", "lods")
  
  for (s in slots_to_rename) {
    if (!is.null(fit2[[s]])) {
      rownames(fit2[[s]]) <- gene_symbols
    }
  }
  
  if (!is.null(fit2$Amean)) {
    names(fit2$Amean) <- gene_symbols
  }
  
  # Also store them in fit2$genes (keep IDs as rownames, symbols as a column)
  fit2$genes <- data.frame(
    ID     = gene_ids,
    Symbol = gene_symbols,
    stringsAsFactors = FALSE,
    row.names = gene_ids  # original IDs are unique, so no duplicate row.names error
  )
  
  #########################
  # 3b. Extract dataframes *no TMM*
  #########################
  
  lfc_matrix <- fit2$coefficients
  
  pval_matrix <- fit2$p.value 
  
  bh_matrix <- apply(pval_matrix, 2, p.adjust, method = "BH")
  rownames(bh_matrix) <- rownames(pval_matrix)
  
  lfc_bh_matrix <- -log10(bh_matrix) * sign(lfc_matrix)
  rownames(lfc_bh_matrix) <- rownames(lfc_matrix)
  
  #########################
  # 4b. Transform and load matrices as .csvs and write *no TMM*
  #########################
  
  print('Write lfc matrix...')
  write.csv(
    t(lfc_matrix),
    file.path(script_dir, "..", data_dir,
              "raw_counts_filtered_logCPM_limma_dmso_lfc_symbols.csv")
  )
  
  print('Write pval matrix...')
  write.csv(
    t(pval_matrix),
    file.path(script_dir, "..", data_dir,
              "raw_counts_filtered_logCPM_limma_dmso_pval_symbols.csv")
  )
  
  print('Write bh matrix...')
  write.csv(
    t(bh_matrix),
    file.path(script_dir, "..", data_dir,
              "raw_counts_filtered_logCPM_limma_dmso_bh_symbols.csv")
  )

} # end if none / all

