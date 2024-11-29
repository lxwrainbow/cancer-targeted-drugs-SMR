file_path1 <- "../MR/OUT/UKB/***"
data1 <- read_delim(file = gzfile(file_path1), delim = "\t")

merged_data <- merge(data, data1, by = "variant")


data2 <- merged_data %>% dplyr::select("rsid","alt","ref","minor_AF.y","beta","se","pval","n_complete_samples")

colnames(data2) <- c("SNP","A1","A2","freq","b", "se","p","n")


data3 <- na.omit(data2)
rm(data1,data2,merged_data)


new_data <- data3 %>% separate_rows(SNP, sep = ",")


write.table(data3, file = "***.txt", sep = "\t", row.names = FALSE, quote = FALSE)
