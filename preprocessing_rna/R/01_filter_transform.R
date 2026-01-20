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
parser <- ArgumentParser(description = 'Filter and transform count data')
parser$add_argument(
  "--data-dir", type = "character", default = "data",
  help = 'Directory name for reading/writing data files (default: "data")'
)
args <- parser$parse_args()

#########################
# 1. Extract data
#########################

script_dir  <- dirname(get_script_path())
data_dir    <- args$data_dir

data_df <- as.data.frame(fread(
  file.path(script_dir, "..", data_dir, "raw_counts.csv")
))

meta_df <- as.data.frame(fread(
  file.path(script_dir, "..", data_dir, "metadata.csv"), header = TRUE
))

# Set rownames to gene_ids for data_df (assuming V1 is gene_ids)
rownames(data_df) <- data_df$V1
data_df$V1 <- NULL

# As type, int
counts_df <- round(data_df)            # just in case
counts_df <- as.matrix(counts_df)
storage.mode(counts_df) <- "integer"

### Format meta_df for edgeR / limma ###

# Set rownames to sample_id
rownames(meta_df) <- meta_df$V1

#########################
# 2. Sanity checks
#########################

## Dimension + alignment checks
stopifnot(ncol(counts_df) == nrow(meta_df))
stopifnot(identical(colnames(counts_df), rownames(meta_df)))

## NA checks
if (sum(is.na(counts_df)) > 0) {
  warning("counts_df contains NA values.")
} else {
  cat("counts_df: no NA values detected.\n")
}

if (sum(is.na(meta_df)) > 0) {
  warning("meta_df contains NA values.")
} else {
  cat("meta_df: no NA values detected.\n")
}

## Check counts look like counts
if (any(counts_df < 0, na.rm = TRUE)) {
  stop("Negative values found in counts_df.")
} else {
  cat("counts_df: no negative values.\n")
}

non_integer <- any(abs(counts_df - round(counts_df)) > 1e-6, na.rm = TRUE)

if (non_integer) {
  warning("Some values in counts_df are not integers.")
} else {
  cat("counts_df: all values appear integer-like.\n")
}

## Library size check
lib_sizes <- colSums(counts_df)
print(summary(lib_sizes))

## Duplicated IDs checks
if (any(duplicated(colnames(counts_df)))) {
  warning("Duplicated sample IDs found in counts_df columns.")
} else {
  cat("No duplicated sample IDs found.\n")
}

if (any(duplicated(rownames(counts_df)))) {
  warning("Duplicated gene IDs found in counts_df rows.")
} else {
  cat("No duplicated gene IDs found.\n")
}

## Drop any groups with < 3 replicates
rep_counts <- table(meta_df$Compound_Concentration)
low_rep_groups <- names(rep_counts)[rep_counts < 3]

if (length(low_rep_groups) > 0) {
  warning(paste(
    "Dropping groups with <3 replicates:", 
    paste(low_rep_groups, collapse = ", ")
  ))
  
  keep_samples <- !(meta_df$Compound_Concentration %in% low_rep_groups)
  meta_df   <- meta_df[keep_samples, , drop = FALSE]
  counts_df <- counts_df[, keep_samples, drop = FALSE]
  
} else {
  cat("All groups have >= 3 replicates. No samples dropped.\n")
}

#########################
# 3. Filter low expressors w/ edgeR (DGEList)
#########################

meta_df$cell_line <- gsub("-", "_", meta_df$cell_line)

# Create condition col (cell_line, drug, dose, time)
meta_df$condition <- paste(
  meta_df$cell_line, meta_df$drug, meta_df$dose, meta_df$time, sep = "_"
)

# Exlore replicates
print(table(meta_df$condition))
cat(table(meta_df$condition))

### CREATE GROUP FACTOR FOR CONTRASTS ###
group <- factor(meta_df$condition)
#group <- relevel(group, ref = "DMSO_0")
levels(group)

y <- DGEList(counts = counts_df)

################################################
#### Filter by minimum count  ####
################################################
keep <- filterByExpr(y, 
                     group = group,
                     min.prop = 1,
                     min.count = 15)

# View genes by min count filter
table(keep)

y <- y[keep, , keep.lib.sizes = FALSE]

## Library size check
lib_sizes <- colSums(y$counts)
print(summary(lib_sizes))

# Write filtered *raw* count matrix (shared for all workflows)
write.csv(
  t(y$counts),
  file.path(script_dir, "..", data_dir, "raw_counts_filtered.csv")
)

#########################
# 4. TMM normalization + Log2-CPM Transformation 
#########################

cat("Creating TMM + Log2CPM transformation...\n")

y_tmm <- y
y_tmm <- calcNormFactors(y_tmm, method = "TMM")

logCPM_tmm <- cpm(y_tmm, log=TRUE, prior.count=1)

write.csv(
  t(logCPM_tmm),
  file.path(script_dir, "..", data_dir,
            "raw_counts_filtered_TMM_logCPM.csv")
)


#########################
# 5. Write meta_df with group col for limma design 
#########################

write.csv(
  meta_df,
  file.path(script_dir, '..', data_dir,
            'metadata_limma_design.csv'),
)

