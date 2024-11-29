setwd("../MR/exp/eQTL/smr_win/smr-1.3.1-win-x86_64")
data <- read.table("chr*_pleio.txt", header = TRUE, sep = "\t")
target_rows <- subset(data, SNP == "top cis_eQTL" & Probe_bp < BP+1000000 & Probe_bp > BP-1000000)
write.table(target_rows, "top cis_eQTL.txt", sep = "\t", quote = FALSE, row.names = FALSE)


setwd("../MR/exp/eQTL/smr_win/smr-1.3.1-win-x86_64")
data <- read.table("chr*.txt", header = TRUE, sep = "\t")
target_rows <- subset(data, Gene %in% c("GeneID*","GeneID*"))
write.table(target_rows, "Gene_pleio.txt", sep = "\t", quote = FALSE, row.names = FALSE)