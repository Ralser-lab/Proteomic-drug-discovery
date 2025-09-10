# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

install.packages("data.table")
install.packages("iq")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")

library(iq)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  setwd(args[1])
  file_path <- args[2]
  out_path <- args[3]
  convert <- as.logical(args[4])
} else {
  stop("Not enough arguments provided")
}

pasef <- fread(file_path, head = TRUE, sep = '\t', dec = '.')

colnames(pasef) <- c('id', 'sample_list', 'protein_list', 'quant')

maxlfq = iq::fast_MaxLFQ(pasef)

pasef_lfq = as.data.frame(t(maxlfq$estimate))

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

install.packages("data.table")
install.packages("iq")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)

library(iq)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  setwd(args[1])
  file_path <- args[2]
  out_path <- args[3]
  convert <- as.logical(args[4])
} else {
  stop("Not enough arguments provided")
}

pasef <- fread(file_path, head = TRUE, sep = '\t', dec = '.')

colnames(pasef) <- c('id', 'sample_list', 'protein_list', 'quant')

maxlfq = iq::fast_MaxLFQ(pasef)

pasef_lfq = as.data.frame(t(maxlfq$estimate))

if (convert) {
  # Extract UniProt IDs from column names
  uniprot_ids <- colnames(pasef_lfq)
  
  # Use AnnotationDbi and org.Hs.eg.db to map UniProt IDs to gene symbols
  gene_symbols <- mapIds(org.Hs.eg.db, keys = uniprot_ids, column = "SYMBOL", keytype = "UNIPROT", multiVals = function(x) paste(x, collapse = ";"))
  
  # Rename the columns in the data frame using the mapping
  colnames(pasef_lfq) <- ifelse(is.na(gene_symbols), colnames(pasef_lfq), gene_symbols)
}

write.csv2(pasef_lfq, out_path)

write.csv2(pasef_lfq, out_path)
