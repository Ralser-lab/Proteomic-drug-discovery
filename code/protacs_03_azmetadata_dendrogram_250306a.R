# Define required CRAN and Bioconductor packages
cran_packages <- c("ape", "ggplot2", "RColorBrewer", "ggnewscale")
bioc_packages <- c("ggtree", "ggtreeExtra")

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

# set absolute path to /data folder in clone github repo in local dir
abspath = '/Users/shaon/Desktop/PROTACS/github_deposition/data/'
setwd(abspath)

## Plots dendrogram of AZ molecules using hierarchical clusters from cleaned AZ metadata
encoded_data = read.csv2('AZcompound_metadata_onehotencoded_240611a.tsv', sep = ',', row.names = 1, 
                        check.names = FALSE) 

# Split metadata into chemical labels and one hot encoded features

chemvars <- names(encoded_data) %in% c('Cluster', 'Drug_Type', "Ligase", 'Target')

logical_data <- encoded_data[!chemvars]

label_data <- encoded_data[chemvars]

logical_data <- as.data.frame(lapply(logical_data, function(column) column == 'True'))

row.names(logical_data) = row.names(encoded_data)

row.names(label_data) = row.names(encoded_data)


# Convert into matrix form , cluster

logical_matrix <- as.matrix(logical_data)

dist_matrix <- dist(logical_matrix, method = "euclidean")

hc <- hclust(dist_matrix, method = "ward.D2")

phylo_tree <- as.phylo(hc)

# Heatmap label

mapdata <- as.data.frame(sapply(label_data, as.character))

mapdata$AZ <- rownames(label_data)

mapdata$Drug_Type <- factor(mapdata$Drug_Type, levels = c('AR-PROTAC', 'Non-PROTAC', 'Txn-PROTAC'))

mapdata$Ligase <- factor(mapdata$Ligase, levels = c('CRBN_lenalidomide_5N', 'CRBN_thalidomide_5N', 'CRBN_dihydrouracil', 'CRBN_thalidomide_6N', 'VHL_amide_tBu', 'Other', 'None'))

mapdata$Ligase <- factor(mapdata$Ligase, levels = c('CRBN_lenalidomide_5N', 'CRBN_dihydrouracil', 'CRBN_thalidomide_6N','CRBN_thalidomide_5N',  'VHL_amide_tBu', 'Other', 'None'))

mapdata$Target <- factor(mapdata$Target, levels = c('AR_piperidine','AR_indole','Other','None'))

# Plot

figout = paste(abspath, '../figures/', sep = '')

setwd(figout)

pdf("ChemicalSeries_dendrogram_plot.pdf", width = 20, height = 10)

options(repr.plot.width = 20, repr.plot.height =10) 

p <- ggtree(phylo_tree, layout = 'circular', size = 0.5)

w = 0.5

p1 <- p + new_scale_fill() + geom_fruit(mapdata, geom = geom_tile, mapping = aes(fill = Drug_Type, y = AZ),
                                         width = w, pwidth=0.1, offset=0.1, inherit.aes = FALSE) + 
  scale_fill_manual(name="Class", values=rev(brewer.pal(3, "Greys")),
                    guide=guide_legend(keywidth=0.5,
                                       keyheight=0.5,
                                       order=1))

p2 <- p1 + new_scale_fill() + geom_fruit(mapdata, geom = geom_tile, mapping = aes(fill = Target, y = AZ),
                                        width = w, pwidth=0.1, offset=0.1) + 
  scale_fill_manual(name="Ligand", values=rev(brewer.pal(5, "Blues")),
                    guide=guide_legend(keywidth=0.5,
                                       keyheight=0.5,
                                       order=2)) 


p3 <- p2 + new_scale_fill() + geom_fruit(mapdata, geom = geom_tile, mapping = aes(fill = Ligase, y = AZ),
                                        width = w, pwidth=0.1, offset=0.1) + 
  scale_fill_manual(name="Degrader", values=rev(brewer.pal(7, "Purples")),
                    guide=guide_legend(keywidth=0.5,
                                       keyheight=0.5,
                                       order=2)) 


custom_spectral <- c('darkblue','blue3','blue2','blue','lightblue','aquamarine','aquamarine2','darkseagreen2','greenyellow','darkolivegreen1','yellow2','darkgoldenrod1','orange','darkorange','red')


mapdata$Cluster <- factor(mapdata$Cluster, levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))

p4 <- p3 + new_scale_fill() + geom_fruit(mapdata, geom = geom_tile, mapping = aes(fill = Cluster, y = AZ),
                                         width = w, pwidth=0.1, offset=0.1) + 
  scale_fill_manual(name="Cluster", values=custom_spectral,
                    guide=guide_legend(keywidth=4,
                                       keyheight=0.5,
                                       ncol = 2,
                                       order=3)) 

p4 + theme(
  plot.margin = unit(c(0, 0, 0, 0), "cm"),
  legend.background=element_rect(fill=NA),
  legend.title=element_text(size=20), 
  legend.text=element_text(size=14),
  legend.spacing.y = unit(0.02, "cm"))

print(p4)

dev.off()


