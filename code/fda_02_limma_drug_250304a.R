
# Define required packages
required_packages <- c(
  "BiocManager", "limma", "EnhancedVolcano", "factoextra", 
  "ggplot2", "dplyr", "tidyr", "pheatmap", "grid", "cowplot", 
  "stats", "ggfortify", 'renv'
)

# Function to install missing packages
install_if_missing <- function(packages) {
  installed <- rownames(installed.packages()) # Get installed packages
  missing <- packages[!(packages %in% installed)] # Check missing ones
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org/")
  }
}

# Ensure Bioconductor packages are installed correctly
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Install missing CRAN packages
install_if_missing(setdiff(required_packages, c("limma", "EnhancedVolcano")))

# Install missing Bioconductor packages
bioc_packages <- c("limma", "EnhancedVolcano")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Load all libraries
lapply(required_packages, library, character.only = TRUE)

message("All required packages are installed and loaded successfully!")

# setwd to correct relative path
setwd('/Users/shaon/Desktop/PROTACS/github_deposition/data')

metadata <- read.csv2('SB_FDA_metadata_250304a.tsv', sep = ',', row.names = 1)

expressionmatrix <- read.csv2('SB_FDA_prmatrix_filtered_50_imputed_50_ltrfm_batched_summarized_250304.tsv', sep = ',', dec = '.', row.names = 1)

expressionmatrix <- t(expressionmatrix)

rownames(expressionmatrix)[rownames(expressionmatrix)=='CYP3A7.CYP3A7.CYP3A51P']<- 'CYP3A7'

metadata[metadata==''] <- 'DMSO'

colnames(metadata) <- c('Plate_Well_', 'Batch_','Plate_','Well_','Content_','Drug_','Cat_','Cas_','Target_')

# Create design contrasts

design <- model.matrix(~0+Drug_,metadata)

# Make target labels compatible with R row annotation

colnames(design) = make.names(colnames(design))

metadata$Drug_ <- make.names(metadata$Drug_)

# Linear models by design matrix

fit <- lmFit(expressionmatrix, design)

# Linear models by contrast matrix

metadata$Drug_ <- as.factor(metadata$Drug_)

# Drop DMSO self contrast
drugs <- levels(metadata$Drug_)[levels(metadata$Drug_) != "DMSO"]

# Define contrasts
drugs <- paste0("Drug_", drugs)

contrast_formulas <- setNames(sapply(drugs, function(drug) paste0(drug, " - Drug_DMSO"), USE.NAMES = FALSE), 
                              sapply(drugs, function(drug) paste0(drug, "vsDrug_DMSO")))

contrast_formulas <- contrast_formulas[names(contrast_formulas) != "Drug_vsDrug_DMSO"]

contrasts <- makeContrasts(contrasts = as.list(contrast_formulas), levels = design)

fit2 <- contrasts.fit(fit, contrasts)

fit3 <- eBayes(fit2, trend = TRUE)

plotSA(fit2, main="Mean-variance trend")
plotMA(fit2, main="Mean-variance trend")
plotMDS(fit2, main="Mean-variance trend")

pval <- fit3$p.value
adjpval <- apply(pval, 2, function(x) p.adjust(x, method = "BH"))
logFC_matrix2 <- fit3$coefficients

signedpval <- -log10(pval) * sign(logFC_matrix2)
signed_adjpval<- -log10(adjpval) * sign(logFC_matrix2)

write.csv2(signedpval, 'FDA_LimmaMatrix_250304a.csv')
write.csv2(signed_adjpval, 'FDA_adjLimmaMatrix_250304a.csv')

pheatmap(logFC_matrix2, scale = 'row',
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
        color=colorRampPalette(c("navy","blue", "lightblue", "white", "pink", "red", "darkred"))(50),
         annotation_col = NULL,  # If you have additional metadata for contrasts, you can add it here
         show_rownames = FALSE,
          
         show_colnames = FALSE)

plot <- function(x, color, ID, plim_y){
  
  rownames(x)[rownames(x)=='CYP3A7.CYP3A7.CYP3A51P']<- 'CYP3A7'
  
  
  keyvals <- ifelse(
    x$logFC < -0.58 & x$adj.P.Val < 0.05, 'blue',
    ifelse(x$logFC > 0.58 & x$adj.P.Val < 0.05, 'red',
           color))
  keyvals[is.na(keyvals)] <- color
  names(keyvals)[keyvals == 'red'] <- '+ 50%'
  names(keyvals)[keyvals == color] <- ''
  names(keyvals)[keyvals == 'blue'] <- '- 50%'
  
  
  p <- EnhancedVolcano(x,lab = rownames(x),
                       x = 'logFC',
                       y = 'adj.P.Val',
                       xlim = c(-1.5, 1.5),
                       ylim = c(-0,plim_y),
                       selectLab = c(ID),
                       ylab = NULL,
                       xlab = NULL,
                       labSize = 8,
                       title = NA,
                       pCutoff = 0.05,
                       axisLabSize = 20,
                       FCcutoff = 0.584,
                       pointSize = 3.0,
                       drawConnectors = TRUE,
                       widthConnectors = 0.75,
                       caption = NA,
                       gridlines.major = FALSE,
                       gridlines.minor = FALSE, hlineWidth = 0.2, vlineWidth = 0.2,
                       border = "full", borderWidth = 0.2, borderColour = "black",
                       colCustom = keyvals,
                       colAlpha = 0.2)
  return(p) 
}


plot2 <- function(x, color, ID, plim_y){
  p <- plot(x, color, ID, plim_y)
  p <- p + theme(axis.ticks = element_line(size = 0.2))
  p <- p + labs(title = NULL, subtitle = NULL)
  p <- p + theme(plot.title = element_text(size = 20))
  p <- p + labs(caption = NULL)
  return(p)
}

result_table1 <- topTable(fit3, coef="Drug_Methotrexate - Drug_DMSO", number=Inf) # good
result_table2 <- topTable(fit3, coef="Drug_Fulvestrant - Drug_DMSO", number=Inf) # good
result_table3 <- topTable(fit3, coef="Drug_Asenapine.maleate - Drug_DMSO", number=Inf) # good
result_table4 <- topTable(fit3, coef="Drug_Doxorubicin..Adriamycin..HCl - Drug_DMSO", number=Inf) # good
result_table5 <- topTable(fit3, coef='Drug_Epirubicin.HCl - Drug_DMSO', number = Inf)
result_table6 <- topTable(fit3, coef="Drug_Tamoxifen - Drug_DMSO", number=Inf) # good
result_table7 <- topTable(fit3, coef="Drug_Ibuprofen. - Drug_DMSO", number = Inf) 
result_table8 <- topTable(fit3, coef="Drug_Clotrimazole - Drug_DMSO", number = Inf)
result_table9 <- topTable(fit3, coef='Drug_Doripenem.Hydrate - Drug_DMSO', number = Inf)
result_table10 <- topTable(fit3, coef='Drug_Amonafide - Drug_DMSO', number = Inf)
result_table11 <- topTable(fit3, coef='Drug_Fluorouracil..5.Fluoracil..5.FU. - Drug_DMSO', number = Inf)

targets <- c('DHFR', 'TK1', 'CYP3A7', 'COX2')

p1 <- plot2(result_table1, 'grey', targets, 4) + labs(title = 'Methotrexate \n(DHFR)') # good
p2 <- plot2(result_table2, 'grey', targets, 4) + labs(title = 'Fulvestrant') # good
p3 <- plot2(result_table3, 'grey', targets, 4) + labs(title = 'Asenapine')
p4 <- plot2(result_table4, 'grey', targets, 4) + labs(title = 'Doxorubicin \n(TK1)')
p5 <- plot2(result_table5, 'grey', targets, 4) + labs(title = 'Epirubicin')
p6 <- plot2(result_table6, 'grey', targets,  4) + labs(title = 'Tamoxifen \n(ER)')
p7 <- plot2(result_table7, 'grey', targets, 4) + labs(title = 'Ibuprofen \n(COX2)')
p8 <- plot2(result_table8, 'grey', targets, 4) + labs(title = 'Clotrimazole \n(CYP450)')
p9 <- plot2(result_table9, 'grey',targets, 4) + labs(title = 'Doripenem')
p10 <- plot2(result_table10, 'grey', targets,4) + labs(title = 'Amonafide')
p11 <- plot2(result_table11, 'grey', targets,4) + labs(title = 'Fluorouracil')

write.csv2(metadata, 'FDA_LimmaMetadata_250304a.csv')

write.csv2(result_table1, 'Methotrexate_DE_250304a.csv')

setwd('/Users/shaon/Desktop/PROTACS/github_deposition/figures')

png('volcanoes_top3_FDA.png', width = 1000, height = 500)
plot_grid(p7, p1, p4, p8, nrow = 1)
dev.off()



