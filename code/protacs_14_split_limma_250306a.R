# Define required CRAN and Bioconductor packages
cran_packages <- c("ggplot2", "dplyr", "tidyr", "pheatmap", "grid", "cowplot")
bioc_packages <- c("limma", "EnhancedVolcano")

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


#set absolute path to /data folder in github repository 
filepath = '/Users/shaon/Desktop/PROTACS/github_deposition/data/'
setwd(filepath)

metadata <- read.csv2('SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_metaforlimma_240611a.tsv', sep = ',', dec = '.', row.names = 1)

expressionmatrix <- read.csv2('SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv', 
                              sep = ',', dec = '.', row.names = 1, header = 1,
                              check.names = FALSE)

colnames(metadata) <- c('Compound_', 'Concentration_', 'Plate_', 'Check_', 'SpRep_', 'SpBatch_', 'SpPosition_', 'Annotation_', "IC50_Glu_", "IC50_Gal_",
                        'E3_ligase_', 'Target_', 'Drug_Type_', 
                        'Binned_Ligase_', 'Binned_Target_', 'Dend_', 'Cluster_', 'Cluster1_', 'Cluster2_')

# subset drugs in the split for GBDT training

split <- read.csv2('AZcompound_metadata_clustered_split_240611a.tsv', sep = ',')

split_filt <- split[split$GBDT_split_idx=='1.0',]

unique_compounds <- unique(split_filt$to_explode)

unique_compounds <- gsub("-", "_", unique_compounds)

metadata[is.na(metadata)] <- 0

metadata[] <- lapply(metadata, function(col) {
  if (is.factor(col)) {
    as.factor(gsub("-", "_", as.character(col)))
  } else if (is.character(col)) {
    gsub("-", "_", col)
  } else {
    col
  }
})

metadata_10 <- subset(metadata, metadata$Concentration_ == 10)

# Filter metadata_10 to keep rows with Compound_ values present in split
metadata_10_filtered <- metadata_10[metadata_10$Compound_ %in% unique_compounds, ]

metadata_DMSO <- subset(metadata, metadata$Concentration_ == 0)

metadata_10_DMSO <- rbind(metadata_10_filtered, metadata_DMSO)

expression_10_DMSO <- expressionmatrix[, colnames(expressionmatrix) %in% rownames(metadata_10_DMSO)]

# Get the order of column names in expression_10_DMSO
column_order <- colnames(expression_10_DMSO)

# Reorder the rows of metadata_10_DMSO based on the column order
metadata_10_DMSO <- metadata_10_DMSO[column_order, , drop = FALSE]

# Ensure the row names match the column names of expression_10_DMSO
rownames(metadata_10_DMSO) <- column_order

metadata_10_DMSO$Cluster2_ <- as.factor(metadata_10_DMSO$Cluster2_)

metadata_10_DMSO$Plate_ <- as.factor(metadata_10_DMSO$Plate_)

metadata_10_DMSO$Concentration_ <- as.numeric(metadata_10_DMSO$Concentration_)

design <- model.matrix(~0+Cluster2_,metadata_10_DMSO)

fit <- lmFit(expression_10_DMSO, design)

drugs <- levels(metadata_10_DMSO$Cluster2_)[levels(metadata_10_DMSO$Cluster2_) != "DMSO"]
drugs <- paste0("Cluster2_", drugs)

contrast_formulas <- setNames(sapply(drugs, function(drug) paste0(drug, " - Cluster2_DMSO"), USE.NAMES = FALSE), 
                              sapply(drugs, function(drug) paste0(drug, "vsCluster2_DMSO")))

contrast_formulas <- contrast_formulas[names(contrast_formulas) != "Cluster2_vsCluster2_DMSO"]

contrasts <- makeContrasts(contrasts = as.list(contrast_formulas), levels = design)

fit2 <- contrasts.fit(fit, contrasts)
fit3 <- eBayes(fit2, trend = TRUE)

plotSA(fit3, main="Mean-variance trend")
plotMA(fit3, main="Mean-variance trend")
plotMDS(fit3, main="Mean-variance trend")

pval <- fit3$p.value
logFC_matrix2 <- fit3$coefficients
signedpval <- -log10(pval) * sign(logFC_matrix2)

colnames(logFC_matrix2) <- rev(c('AR PROTAC 5N thal', 'AR PROTAC Dihydro U', 'AR PROTAC 5N lena', 'Other AR PROTAC', 'AR PROTAC 6N', 'AR PROTAC VHL amide', 'Non PROTAC', 'TXN PROTAC'))

colnames(pval) <- gsub("Cluster2_", "", colnames(pval))
colnames(pval) <- gsub("- DMSO", "", colnames(pval))
colnames(logFC_matrix2) <- gsub("Cluster2_", "", colnames(logFC_matrix2))
colnames(logFC_matrix2) <- gsub("- DMSO", "", colnames(logFC_matrix2))
colnames(signedpval) <- gsub("Cluster2_", "", colnames(logFC_matrix2))
colnames(signedpval) <- gsub("- DMSO", "", colnames(logFC_matrix2))

result_table1 <- topTable(fit3, coef="Cluster2_1.0 - Cluster2_DMSO", number=Inf) # good
result_table2 <- topTable(fit3, coef="Cluster2_2.0 - Cluster2_DMSO", number=Inf) # good
result_table3 <- topTable(fit3, coef="Cluster2_3.0 - Cluster2_DMSO", number=Inf) # good
result_table4 <- topTable(fit3, coef="Cluster2_4.0 - Cluster2_DMSO", number=Inf) # good
result_table5 <- topTable(fit3, coef="Cluster2_5.0 - Cluster2_DMSO", number=Inf) # good
result_table6 <- topTable(fit3, coef="Cluster2_6.0 - Cluster2_DMSO", number=Inf)
result_table7 <- topTable(fit3, coef="Cluster2_7.0 - Cluster2_DMSO", number=Inf)
result_table8 <- topTable(fit3, coef="Cluster2_8.0 - Cluster2_DMSO", number=Inf)


write.csv2(result_table8, 'Cluster_T5N_Limma_250306a.csv', row.names = TRUE)
write.csv2(result_table7, 'Cluster_URA_Limma_250306a.csv', row.names = TRUE)
write.csv2(result_table6, 'Cluster_L5N_Limma_250306a.csv', row.names = TRUE)
write.csv2(result_table5, 'Cluster_OTH_Limma_250306a.csv', row.names = TRUE)
write.csv2(result_table4, 'Cluster_T6N_Limma_250306a.csv', row.names = TRUE)
write.csv2(result_table3, 'Cluster_VHL_Limma_250306a.csv', row.names = TRUE)
write.csv2(result_table2, 'Cluster_NON_Limma_250306a.csv', row.names = TRUE)
write.csv2(result_table1, 'Cluster_TXN_Limma_250306a.csv', row.names = TRUE)

