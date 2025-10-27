# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

install.packages("data.table")
install.packages('iq')

library(iq)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  setwd(args[1])
  file_path <- args[2]
  out_path <- args[3]
} else {
  stop("Not enough arguments provided")
}

pasef <- fread(file_path, head = TRUE,
                 sep = '\t', dec = '.')

colnames(pasef) <- c('id','sample_list', 'protein_list','quant')

maxlfq = iq::fast_MaxLFQ(pasef)

pasef_lfq = as.data.frame(t(maxlfq$estimate))

write.csv2(pasef_lfq, out_path)
