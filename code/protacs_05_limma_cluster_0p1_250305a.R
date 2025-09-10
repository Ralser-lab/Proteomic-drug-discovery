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


# Performs Limma to generate DE matrices by Cluster ID in PROTAC dataset 0.1 micromolar

# Performs Limma to generate DE matrices by Cluster ID in PROTAC dataset at 10 micromolar
filepath = '/Users/shaon/Desktop/PROTACS/github_deposition/data/'
setwd(filepath)

metadata <- read.csv2('SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_metaforlimma_240611a.tsv', sep = ',', dec = '.', row.names = 1)

expressionmatrix <- read.csv2('SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv', 
                              sep = ',', dec = '.', header = 1, 
                              check.names = FALSE)

test <- subset(expressionmatrix, expressionmatrix[,1]!='')

row.names(test) <- test[,1]

expressionmatrix <- test

test[,1] <- NULL

expressionmatrix <- test

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

metadata_10_DMSO <- subset(metadata, metadata$Concentration_ == 0.1 | metadata$Concentration_ == 0)

#metadata_10_DMSO <- subset(metadata, metadata$Concentration_ == 1 | metadata$Concentration_ == 0)

#metadata_10_DMSO <- subset(metadata, metadata$Concentration_ == 10 | metadata$Concentration_ == 0)

expression_10_DMSO <- expressionmatrix[, colnames(expressionmatrix) %in% rownames(metadata_10_DMSO)]

z <- t(table(metadata_10_DMSO$Compound_, metadata_10_DMSO$Cluster2_))
z

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

adj_pval <- apply(pval, 2, function(x) p.adjust(x, method = "BH"))

logFC_matrix1 <- fit3$coefficients

logFC_matrix2 <- -log10(pval) * sign(logFC_matrix1)

logFC_matrix3 <- -log10(adj_pval) * sign(logFC_matrix1)

write.csv2(logFC_matrix1, 'Cluster_LFC_0p1uM_250305a.csv')

write.csv2(logFC_matrix2, 'Cluster_LFCxPval_0p1uM_250305a.csv')

write.csv2(logFC_matrix3, 'Cluster_LFCxadjPval_0p1uM_250305a.csv')

colnames(pval) <- gsub("Cluster2_", "", colnames(pval))
colnames(pval) <- gsub("- DMSO", "", colnames(pval))
colnames(logFC_matrix2) <- gsub("Cluster2_", "", colnames(logFC_matrix2))
colnames(logFC_matrix2) <- gsub("- DMSO", "", colnames(logFC_matrix2))

plot <- function(x, color){
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
                       xlim = c(-1.5, 1.5),
                       ylim = c(-0,30),
                       ylab = bquote(~-Log[10]~adjusted~italic(P)),
                       title = NA,
                       pCutoff = 0.05,
                       axisLabSize = 20,
                       FCcutoff = 0.58,
                       pointSize = 3.0,
                       caption = NA,
                       legendPosition = "none",
                       gridlines.major = FALSE,
                       gridlines.minor = FALSE, hlineWidth = 0.2, vlineWidth = 0.2,
                       border = "full", borderWidth = 0.2, borderColour = "black",
                       col = c("grey", "grey", "grey", color))
  
  return(p) 
}


plot2 <- function(x, color){
  p <- plot(x, color)
  p <- p + theme(axis.ticks = element_line(size = 0.2))
  p <- p + labs(title = NULL, subtitle = NULL)
  p <- p + theme(plot.title = element_text(size = 18))
  p <- p + labs(caption = NULL)
  return(p)
}


  
result_table1 <- topTable(fit3, coef="Cluster2_1.0 - Cluster2_DMSO", number=Inf) # good
result_table2 <- topTable(fit3, coef="Cluster2_2.0 - Cluster2_DMSO", number=Inf) # good
result_table3 <- topTable(fit3, coef="Cluster2_3.0 - Cluster2_DMSO", number=Inf) # good
result_table4 <- topTable(fit3, coef="Cluster2_4.0 - Cluster2_DMSO", number=Inf) # good
result_table5 <- topTable(fit3, coef="Cluster2_5.0 - Cluster2_DMSO", number=Inf) # good
result_table6 <- topTable(fit3, coef="Cluster2_6.0 - Cluster2_DMSO", number=Inf)
result_table7 <- topTable(fit3, coef="Cluster2_7.0 - Cluster2_DMSO", number=Inf)
result_table8 <- topTable(fit3, coef="Cluster2_8.0 - Cluster2_DMSO", number=Inf)
result_table9 <- topTable(fit3, coef="Cluster2_9.0 - Cluster2_DMSO", number=Inf)
result_table10 <- topTable(fit3, coef="Cluster2_10.0 - Cluster2_DMSO", number=Inf)
result_table11 <- topTable(fit3, coef="Cluster2_11.0 - Cluster2_DMSO", number=Inf)
result_table12 <- topTable(fit3, coef="Cluster2_12.0 - Cluster2_DMSO", number=Inf)
result_table13 <- topTable(fit3, coef="Cluster2_13.0 - Cluster2_DMSO", number=Inf)
result_table14 <- topTable(fit3, coef="Cluster2_14.0 - Cluster2_DMSO", number=Inf)
result_table15 <- topTable(fit3, coef="Cluster2_15.0 - Cluster2_DMSO", number=Inf)


p1 <- plot2(result_table1, 'darkblue') + labs(title = 'Cluster 1, AR-PROTAC: \nVHL-Other') # good
p2 <- plot2(result_table2, 'blue3') + labs(title = 'Cluster 2, AR-PROTAC: \nThalidomide 6N-indole') # good
p3 <- plot2(result_table3, 'blue2') + labs(title = 'Cluster 3, AR-PROTAC: \nOther-Other')
p4 <- plot2(result_table4, 'blue') + labs(title = 'Cluster 4, TXN-PROTAC: \nVHL/Thalidomide 6N-Other')
p5 <- plot2(result_table5, 'lightblue') + labs(title = 'Cluster 5, AR-PROTAC: \nVHL-Indole')
p6 <- plot2(result_table6, 'aquamarine') + labs(title = 'Cluster 6, AR-PROTAC: \nThalidomide 5N-Other')
p7 <- plot2(result_table7, 'aquamarine2') + labs(title = 'Cluster 7, TXN-PROTAC: \nOther-Other')
p8 <- plot2(result_table8, 'darkseagreen2') + labs(title = 'Cluster 8, AR-PROTAC: \nVHL/Thalidomide 6N-Piperidine')
p9 <- plot2(result_table9, 'greenyellow') + labs(title = 'Cluster 9, AR-PROTAC: \nThalidomide 6N-Other')
p10 <- plot2(result_table10, 'darkolivegreen1') + labs(title = 'Cluster 10, AR-PROTAC: \nThalidomide 5N-Piperidine')
p11 <- plot2(result_table11, 'yellow2') + labs(title = 'Cluster 11, AR-PROTAC: \nDihydrouracyl-Piperidine')
p12 <- plot2(result_table12, 'darkgoldenrod1') + labs(title = 'Cluster 12, AR-PROTAC: \nOther-Piperidine')
p13 <- plot2(result_table13, 'orange') + labs(title = 'Cluster 13, Not PROTAC')
p14 <- plot2(result_table14, 'orangered') + labs(title = 'Cluster 14, AR-PROTAC: \nLenalinomide 5N-Other')
p15 <- plot2(result_table15, 'firebrick') + labs(title = 'Cluster 15, AR-PROTAC: \nLenalinomide 5N-Piperidine')

write.csv2(result_table1, 'Cluster1_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table2, 'Cluster2_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table3, 'Cluster3_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table4, 'Cluster4_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table5, 'Cluster5_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table6, 'Cluster6_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table7, 'Cluster7_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table8, 'Cluster8_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table9, 'Cluster9_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table10, 'Cluster10_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table11, 'Cluster11_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table12, 'Cluster12_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table13, 'Cluster13_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table14, 'Cluster14_0p1uM_Limma_250305a.csv', row.names = TRUE)
write.csv2(result_table15, 'Cluster15_0p1uM_Limma_250305a.csv', row.names = TRUE)

setwd(paste(filepath, '../figures', sep=""))
library('cowplot')

png('volcanoes_0p1uM.png', width = 1200, height = 480)
plot_grid(p1, p6, p9, p14, nrow = 1)
dev.off()

