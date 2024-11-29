
library(TwoSampleMR)
setwd("../MR/TSMR")

a <- read.table("stroke.txt",header = T)
head(a)
a <- subset(a, p<5e-08)

exp_dat1 <- a %>% dplyr::select("SNP","p") 
colnames(exp_dat1)[1] <- 'rsid' 
colnames(exp_dat1)[2] <- 'pval'

stroke <- ieugwasr::ld_clump_local(exp_dat1,clump_p=0.05,clump_r2=0.01,clump_kb=10000, 
                                   bfile = "../MR/Decode/1kg.v3/EUR", 
                                   plink_bin="../MR/Decode/plink_win64_20230116/plink")

sameSNP <- intersect(a$SNP,stroke$rsid)

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
                                    beta_col = "b",se_col = "se",eaf_col = "freq",
                                    effect_allele_col = "A1",other_allele_col = "A2",
                                    pval_col = "p",samplesize_col = "samplesize",
                                    chr_col = "chr",pos_col = "pos")

f <- (exp_dat$beta.exposure/exp_dat$se.exposure)^2
exp_dat $ f <-f

library(data.table)
library(tidyverse)

file <- c('****.txt.gz')

pQTL <- file 
data0 <- fread(pQTL,data.table = F,integer64 = "numeric")

colnames(data0)

c <- data0 %>% dplyr::select("rsids","effectAllele","otherAllele","ImpMAF","Beta","SE","Pval","N")

colnames(c) <- c("SNP","effect_allele","other_allele","eaf","beta","se", 'p','samplesize')

d <- merge(exp_dat,c,by.x = "SNP",by.y ="SNP")
head(d)

write.csv(d,file = "outcome_dat.csv")

outcome_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "outcome_dat.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p")
View(VEGFA121_out_dat)
dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat = VEGFA121_out_dat)
write.csv(dat,file = "harmonized_dat.csv")

res <- mr(dat)
write.csv(res,file = "mr_dat.csv")
OR <- generate_odds_ratios(res)
View(OR)
write.csv(OR, file = "mr_odds_dat.csv")

library(MRPRESSO)
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
          SignifThreshold = 0.05)

het <- mr_heterogeneity(dat)
het

pleio <- mr_pleiotropy_test(dat)
pleio

single <- mr_leaveoneout(dat)
mr_leaveoneout_plot(single)

mr_scatter_plot(res,dat)

res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)

mr_funnel_plot(res_single)