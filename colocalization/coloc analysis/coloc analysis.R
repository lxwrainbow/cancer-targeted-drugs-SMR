library(coloc)
library(remotes)
library(dplyr)

 setwd("../smr_win/SMR/coloc/")
 eqtl = data.table::fread("./EXP/coloc_Gene.txt",sep = "\t",header=T)
 gwas = data.table::fread("./CVD/CVD_coloc_Gene.txt",sep = "\t",header=T)

 data <- merge(eqtl,gwas,by="SNP")
 data <- data[!duplicated(data$SNP),]
 
 data <- data %>% filter((A1.x==A1.y&A2.x==A2.y)|(A1.x==A2.y&A2.x==A1.y)) 
 data <- data %>% mutate(BETA.y = ifelse(A1.x==A1.y,BETA.y,-BETA.y))

 data$VAR.x = data$SE.x^2
 data$VAR.y = data$SE.y^2
 data = data[data$VAR.x!=0 & data$VAR.y!=0 ,]
 
gwas = data[,c("BETA.y","VAR.y","SNP","MAF.y","N")]
eqtl = data[,c("BETA.x","VAR.x","SNP")]
colnames(gwas)=c("beta","varbeta","snp","MAF","N")
colnames(eqtl)=c("beta","varbeta","snp")
eqtl = as.list(eqtl)
gwas = as.list(gwas)


 gwas$type = "quant"
 eqtl$type = "cc"
 
  res = coloc.abf(eqtl,gwas,p1=1e-4,p2=1e-4,p12=1e-5)
  res$summary
  
  library(devtools) 
  devtools::install_github("boxiangliu/locuscomparer",force = T)


  gwas = cbind(data$SNP,data$P.y)
    colnames(gwas) = c("rsid","pval")
    eqtl = cbind(data$SNP,data$P.x)
    colnames(eqtl) = c("rsid","pval")
    write.table(gwas,"gwas.tsv",col.names = T,row.names = F,sep="\t",quote = F)
    write.table(eqtl,"eqtl.tsv",col.names = T,row.names = F,sep="\t",quote = F)
    locuscompare("gwas.tsv","eqtl.tsv",title1 = "CVD GWAS",title2 = "Gene eQTL",
                 legend = T,snp = "top cis_eQTL")  


  