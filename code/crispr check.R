library(dplyr)
library(tidyverse)
library(dslabs)
data(heights)
aa %>% 
  ggplot(aes(`ACH-000001`)) +
  # filter(.> 0) %>% 
  geom_histogram(binwidth = 0.01, fill = "blue", col = "black") 

CRISPR_gene_effect %>% 
  ggplot(aes(`KRAS (3845)`)) +
  # filter(.> 0) %>% 
  geom_histogram(binwidth = 0.01, fill = "blue", col = "black") 


CRISPR_gene_effect %>% 
  ggplot(aes(`TP53 (7157)`)) +
  # filter(.> 0) %>% 
  geom_histogram(binwidth = 0.01, fill = "blue", col = "black") 

CRISPR_gene_effect %>% 
  ggplot(aes(`CDKN2A (1029)`)) +
  # filter(.> 0) %>% 
  geom_histogram(binwidth = 0.01, fill = "blue", col = "black") 

CRISPR_gene_effect %>% 
  ggplot(aes(`SMAD4 (4089)`)) +
  # filter(.> 0) %>% 
  geom_histogram(binwidth = 0.01, fill = "blue", col = "black") 


CRISPR_gene_effect=CRISPR_gene_effect[,-1]
aa=t(CRISPR_gene_effect)
aa=as.data.frame(aa)
aa[,"ACH-000017"]


library(gamlss)
# library(gamlss.dist)
# library(gamlss.add)
# 例如，设置type = "realline"将尝试在整个实线上定义的所有已实现分布，
# 而type = "realsplus"仅尝试在实正线上定义的分布。
aa[,"ACH-000017"]#"Sinh-Arcsinh"
aa[,"ACH-000164"][aa[,"ACH-000164"]<0]#Family:  c("SEP4", "skew exponential power type 4") 
-aa[,"ACH-000164"][aa[,"ACH-000164"]<0]#c("SEP4", "skew exponential power type 4") 
-aa[,"ACH-000164"][aa[,"ACH-000164"]<0]# realplus 

fit <- fitDist(-aa[,"ACH-000164"][aa[,"ACH-000164"]<0], k = 2, type = "realAll", trace = FALSE, try.gamlss = TRUE)

fit <- fitDist(-aa[,"ACH-000164"][aa[,"ACH-000164"]<0], k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)

summary(fit)

bb=aa[order(aa$`ACH-000164`,decreasing = F),]
rownames(bb)
min(aa[,"ACH-000060"])
aa["RAN (5901)","ACH-000060"]
# [1] "RAN (5901)"         "HSPE1 (3336)"       "RPL15 (6138)"       "POLR2L (5441)"      "SMU1 (55234)"      
# [6] "RPL4 (6124)"        "RPS29 (6235)"       "RPL18A (6142)"      "RPS6 (6194)"        "TXNL4A (10907)"    
# [11] "SNRPB (6628)"       "PRPF38A (84950)"    "SNRNP200 (23020)"   "SNRPF (6636)"       "SNRPA1 (6627)" 
data = rnorm(1000,mean=5, sd=0.75)
res = fitDist(data, k = 2, type = "realAll", trace = FALSE, try.gamlss = TRUE)
summary(res)

##################
CRISPR_gene_effect
eset=aa
eset=na.omit(eset)
sigScore <- function(eset, methods = "PCA") {
  
  if (methods == "PCA") {
    # message(paste0("Calculating siganture score using PCA function"))
    pc <- prcomp(t(eset),center = F,na.action = na.omit, scale. = FALSE)
    sigs <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(eset)))
  } else {
    # message(paste0("Calculating siganture score using mean of signature genes"))
    sigs <- colMeans(eset)
  }
  return(sigs)
}

tol：数值型变量，标准差低于此数值的主成分将被忽略，与rank共同作用。

rank：设置待分析的主成分个数，标准差排名靠后的主成分将被舍弃，与tol共同作用。
