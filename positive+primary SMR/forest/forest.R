
library(forestploter)
library(data.table)

setwd("../MR/exp/eQTL/smr_win/SMR/forest")

forest <- read.table("Gene.txt", sep = "\t", header = T)
forest$p.value <- format(forest$p.value, scientific = F , digits = 2)
forest$HEIDI.Test <- format(forest$HEIDI.Test, scientific = F, digits = 2)


forest$se<- (log(as.numeric(forest$high95)) - log(as.numeric(forest$OR)))/1.96


forest$high95 <-as.numeric(forest$high95)
forest$low95 <-as.numeric(forest$low95)

forest$`OR(95% CI)` <- ifelse(is.na(forest$se), "",
                          sprintf("%.2f(%.2f-%.2f)", forest$OR,forest$low95,forest$high95))#sprintF返回字符和可变量组合

forest$` `<- paste(rep(" ", 20), collapse = " ")


tm <- forest_theme(base_size = 12,  
               
                   ci_pch = 20,
                   ci_col = "black",  
                   ci_fill = "black",
                   ci_alpha = 1,        
                   ci_lty = 1,          
                   ci_lwd = 1.5,  
                   refline_lwd = 1,   
                   refline_lty = "solid",
                   refline_col = "black",
                   title_just = "left",
                   title_cex = 1)

forest(forest[,c(1,8,9,5,6)],
       est = forest$OR,       
       lower = forest$low95,     
       upper = forest$high95,   
       sizes = 0.4,
       ci_column = 3,   
       ref_line = 1,
       arrow_lab = c("Decreased risk", "Increased risk"),
       xlim = c(0.5,1.5),
       ticks_at = c(0.5,0.8,1,1.2,1.5),
       xlab = "OR per 1 SD decrease in gene expression(95%CI)",
       title = "Gene inhibition",
       theme = tm)
