########################exposure myquery文件###################################

setwd("../smr_win/smr-1.3.1-win-x86_64")
data <- read.table("chr*_coloc.txt", header = TRUE, sep = "\t")

target_rows <- subset(data, Gene == "Gene ID")

write.table(target_rows, "coloc_Gene.txt", sep = "\t", quote = FALSE, row.names = FALSE)

######################outcome coloc#############################################
library(readr)
library(tidyverse)
library(TwoSampleMR)

setwd("../SMR/coloc/out")
file_path <- "../SMR/coloc/out/CVD.gz"
data <- read_delim(file = gzfile(file_path), delim = "\t")

 data1 <- data %>% dplyr::select("rsids","#chrom","pos","alt","ref","af_alt","beta","sebeta","pval")

 colnames(data1) <- c("SNP","CHR","BP","A1","A2","MAF","BETA", "SE","P")
 data2 <- na.omit(data1)
data3 <- data2 %>% filter(CHR==*,BP >= top cis_eQTL BP -1000000,BP <= top cis_eQTL BP +1000000)

new_data <- data3 %>% separate_rows(SNP, sep = ",")

write.table(new_data, file = "CVD_coloc_Gene.txt", sep = "\t", row.names = FALSE, quote = FALSE)
