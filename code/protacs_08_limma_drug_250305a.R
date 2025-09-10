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

# Performs Limma to generate DE matrices by Drug ID in PROTAC dataset

# set absolute path to /data in the cloned github repository
filepath = '/Users/shaon/Desktop/PROTACS/github_deposition/data/'
setwd(filepath)

metadata <- read.csv2('SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_metaforlimma_240611a.tsv', sep = ',', dec = '.', row.names = 1)

expressionmatrix <- read.csv2('SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv', 
                              sep = ',', dec = '.', row.names = 1, header = 1,
                              check.names = FALSE)

colnames(metadata) <- c('Compound_', 'Concentration_', 'Plate_', 'Check_', 'SpRep_', 'SpBatch_', 'SpPosition_', 'Annotation_', "IC50_Glu_", "IC50_Gal_",
                        'E3_ligase_', 'Target_', 'Drug_Type_', 'Binned_Ligase_', 'Binned_Target_', 'Cluster_', 'Dend_','Cluster2_')

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

metadata_10_DMSO <- subset(metadata, metadata$Concentration_ == 0.1 | metadata$Concentration_ == 1 |
                         metadata$Concentration_ == 10 | metadata$Concentration_ == 0)

metadata_10_DMSO <- subset(metadata, metadata$Concentration_ == 10 | metadata$Concentration_ == 0)

expression_10_DMSO <- expressionmatrix[, colnames(expressionmatrix) %in% rownames(metadata_10_DMSO)]

z <- t(table(metadata_10_DMSO$Compound_, metadata_10_DMSO$Cluster2_))
z

metadata_10_DMSO$Cluster2_ <- as.factor(metadata_10_DMSO$Cluster2_)

metadata_10_DMSO$Plate_ <- as.factor(metadata_10_DMSO$Plate_)

indices <- metadata_10_DMSO$Cluster2_ == 'DMSO'  # Find indices where Cluster2_ is 'DMSO'

metadata_10_DMSO$Annotation_[indices] <- 'DMSO' 

metadata_10_DMSO$Compound_[indices] <- 'DMSO' 

metadata_10_DMSO$Annotation_ <- as.factor(metadata_10_DMSO$Annotation_)

metadata_10_DMSO$Compound_ <- as.factor(metadata_10_DMSO$Compound_)

# Format Annotation to be compatible with R dataframe
levels(metadata_10_DMSO$Annotation_) <- make.names(levels(metadata_10_DMSO$Annotation_))

design <- model.matrix(~0+Compound_,metadata_10_DMSO)

fit <- lmFit(expression_10_DMSO, design)

drugs <- levels(metadata_10_DMSO$Compound)[levels(metadata_10_DMSO$Compound) != "DMSO"]
drugs <- paste0("Compound_", drugs)

contrast_formulas <- setNames(sapply(drugs, function(drug) paste0(drug, " - Compound_DMSO"), USE.NAMES = FALSE), 
                              sapply(drugs, function(drug) paste0(drug, "vsCompound_DMSO")))

contrasts <- makeContrasts(contrasts = as.list(contrast_formulas), levels = design)

fit2 <- contrasts.fit(fit, contrasts)

fit3 <- eBayes(fit2, trend = TRUE)

pval <- fit3$p.value

adjpval <- apply(pval, 2, function(x) p.adjust(x, method = "BH"))

logFC_matrix2 <- fit3$coefficient

limma_matrix2 <- -log10(pval) * sign(logFC_matrix2)

limma_matrix3 <- -log10(adjpval) * sign(logFC_matrix2)

#write.csv2(logFC_matrix2, 'Drug_LFC_250305a.csv')

write.csv2(limma_matrix2, 'Drug_LFCxPval_250305a.csv')

write.csv2(limma_matrix2, 'Drug_LFCxadjPval_250305a.csv')

plot <- function(x, color, vector, y_lim){
  ordered <- x[order(x$logFC,decreasing=TRUE),]
  ordered$rank <- 1:nrow(ordered)
  
  for(i in 1:nrow(ordered)){
    if(ordered$logFC[i] > 0){
      ordered$color[i] <- 'PINK'
    }
    else{
      ordered$color[i] <- 'LIGHTBLUE'
    }
  }
  
  for(i in 1:nrow(ordered)) {
    if(i < 10) {
      ordered$names[i] <- rownames(ordered)[i]
      ordered$color[i] <- 'RED'
    } else if (i > nrow(ordered)-11){
      ordered$names[i] <- rownames(ordered)[i]
      ordered$color[i] <- 'BLUE'
    }
    else {
      ordered$names[i] <- NA
    }
  }
  
  p <- ggplot(ordered, aes(x = rank, y = logFC, label = names, color = color)) +
    geom_point() + geom_text_repel() + theme_bw() + ylim(-5,5) +  scale_color_identity() +
    geom_hline(yintercept = 0, linetype = 'dashed')
  
  p <- EnhancedVolcano(x,lab = rownames(x),
                       x = 'logFC',
                       y = 'adj.P.Val',
                       xlim = c(-3.0, 3.0),
                       ylim = c(-0,y_lim),
                       title = NA,
                       pCutoff = 0.05,
                       axisLabSize = 20,
                       FCcutoff = NA,
                       pointSize = 3.0,
                       boxedLabels = TRUE,
                       drawConnectors = TRUE,
                       caption = NA,
                       legendPosition = "none",
                       gridlines.major = FALSE,
                       gridlines.minor = FALSE, hlineWidth = 0.2, vlineWidth = 0.2,
                       border = "full", borderWidth = 0.2, borderColour = "black",
                       col = c("grey", "grey", color, color))
  
  return(p) 
}

plot2 <- function(x, color, vector, y_lim){
  p <- plot(x, color, vector, y_lim)
  p <- p + theme(axis.ticks = element_line(size = 0.2))
  p <- p + labs(title = NULL, subtitle = NULL)
  p <- p + theme(plot.title = element_text(size = 14))
  p <- p + labs(caption = NULL)
  return(p)
}


