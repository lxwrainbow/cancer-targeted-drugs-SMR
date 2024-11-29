
library(TwoSampleMR)

setwd("../MR/TSMR")


a <- read.table("Gene.txt",header = T)

exp_dat1 <- a %>% dplyr::select("SNP","p") 

colnames(exp_dat1)[1] <- 'rsid' 
colnames(exp_dat1)[2] <- 'pval'

FGFR1 <- ieugwasr::ld_clump_local(exp_dat1,clump_p=0.05,clump_r2=0.01,clump_kb=10000, 
                                  bfile = "../MR/Decode/1kg.v3/EUR", 
                                  plink_bin="../MR/Decode/plink_win64_20230116/plink")

sameSNP <- intersect(a$SNP,Gene$rsid)

if (length(sameSNP) == 0) {
  print(paste("Warning: No intersecting SNPs found in file", file))
  next 
}
data <- a[a$SNP %in% sameSNP,] %>% dplyr::arrange(p) %>% dplyr::distinct(SNP,.keep_all = T)


if (nrow(data) == 0) {
  print(paste("Warning: No data remaining after filtering with intersecting SNPs in file", file))
  next
}

exp_dat <- TwoSampleMR::format_data(data, type = "exposure",snp_col  = "SNP",
                                    beta_col = "b",se_col = "SE",eaf_col = "Freq",
                                    effect_allele_col = "A1",other_allele_col = "A2",
                                    pval_col = "p")

f <- (exp_dat$beta.exposure/exp_dat$se.exposure)^2
exp_dat $ f <-f


library(MRPRESSO)
library(data.table)
library(tidyverse)


c <- read.table("CVD.txt", header=T)
View(c)
d <- merge(exp_dat,c,by.x = "SNP",by.y ="SNP")
View(d)
write.csv(d,file = "CVD_out_dat.csv")
CVD_out_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "CVD_out_dat.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p")

dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat = CVD_out_dat)
write.csv(dat,file = "harmonized_dat.csv")

res <- mr(dat)
write.csv(res,file = "mr_dat.csv")
OR <- generate_odds_ratios(res)
View(OR)
write.csv(OR, file = "mr_odds_dat.csv")




