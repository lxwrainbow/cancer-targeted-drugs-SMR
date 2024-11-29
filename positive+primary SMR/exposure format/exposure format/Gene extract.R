
data <- read.table("chr*.txt", header = TRUE, sep = "\t")
target_rows <- subset(data, Gene == "GeneID")
target_rows$Gene <- "GeneName"
write.table(target_rows, "Gene.txt", sep = "\t", quote = FALSE, row.names = FALSE)