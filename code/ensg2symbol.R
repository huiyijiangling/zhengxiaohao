dataset=ICGC_PACA_CA_seq_raw
dataset[is.na(dataset)] <- 0
load("C:/Users/zxh/Desktop/R/gtf/exonlength.Rdata")


gen <- intersect(rownames(dataset),rownames(exons_gene_len))
dataset_tpm=dataset[gen,]
eff_length=exons_gene_len[gen,]
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

dataset_tpm <- apply(dataset_tpm,2,function(x) countToTpm(x,eff_length$eff_length))
colSums(dataset_tpm)
###################################### 4check
median(dataset_tpm["ENSG00000143947",])
median(rnatpm_Pancreas["ENSG00000143947",])
log2(dataset_tpm["ENSG00000143947",]+0.001)
rnatpmlog2001_Pancreas["ENSG00000143947",]

dataset=dataset_tpm


library(AnnoProbe)
probe2gene=data.frame(probe_id=rownames(dataset),symbol=rownames(dataset))
source("./updateName.R")
forthetree=ensg2SymbolUsingNCBImodify(alias=probe2gene$symbol,gene.info.file = "Homo_sapiens.gene_info.gz")
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
if(F){
  forthetree=updateName(probe2gene$symbol)
  table(forthetree$gene_biotype)
  probe2gene=merge(forthetree,probe2gene,by.x="ALIAS",by.y="symbol",all.y=T)
}
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
probe2gene=unique(probe2gene[,c("probe_id","Symbol")])
# probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
colnames(probe2gene)=c("probe_id","symbol")
probes_expr=dataset
genes_expr <- filterEM(probes_expr,probe2gene)
ICGC_PACA_CA_seq_raw_symbol_tpm=genes_expr
write.csv(ICGC_PACA_CA_seq_raw_symbol_tpm,"ICGC_PACA_CA_seq_raw_symbol_tpm.csv",quote = T)
read.csv("ICGC_PACA_CA_seq_raw_symbol_tpm.csv")


#载入R包；
library(ggcorrplot)
library(ggplot2)
#直接快速绘制整个相关性热图；

# quickcor(mtcars, cluster = TRUE,cor.test = TRUE) +
#   geom_colour() +
#   geom_mark(size=3,color="white",fontface=1)+
#   scale_fill_gradientn(colours = c("#77C034","white" ,"#C388FE"))+
#   geom_panel_grid(colour = "white",size = 1)
icgc_CA_tpm_immune=read.csv("C:/Users/zxh/Desktop/R/ftoimmune/xCell_ICGC_PACA_CA_seq_raw_symbol_tpm_xCell_2209060521.csv")

rownames(icgc_CA_tpm_immune)=icgc_CA_tpm_immune[,1]
icgc_CA_tpm_immune=icgc_CA_tpm_immune[,-1]
ICGC_PACA_CA_seq_raw_symbol_tpm["FTO",]
cor=rbind(ICGC_PACA_CA_seq_raw_symbol_tpm["FTO",],icgc_CA_tpm_immune)

library(tidyverse)
library(corrplot)
# p.mat=res1$p,
data(mtcars)
cor=t(cor)
cor=as.data.frame(cor)
# cor=as.numercoric(cor)
res1 <- cor.mtest(cor, conf.level = .95, method = "spearman")
cor(cor) %>%   corrplot(method = "square",
                        type = "lower",tl.srt = 45,tl.col = "black")

pp1<- cor(cor, method ="spearman") %>%  corrplot(type="lower",method = "square",
                       insig="label_sig",p.mat = res1$p,outline = "white",
                       sig.level=c(.001,.01,.05),
                       pch.cex=.9,pch.col="green",tl.srt=45,tl.col="black")
pp2=res1[["p"]]
pp2=as.data.frame(pp2)
pp1=as.data.frame(pp1)
pp1
ppp_icgc_ca=cbind(pp1[,1],pp2[,1])
rownames(ppp_icgc_ca)=c("FTO",rownames(icgc_CA_tpm_immune))
ppp_icgc_ca_lose3=ppp_icgc_ca[-c(66,67,68),]
write.csv(ppp_icgc_ca_lose3,file="ppp_icgc_ca_lose3.csv")
save(ppp_icgc_ca_lose3,ppp_icgc_ca,file="ppp_icgc_ca_lose3.Rdata")

load("ICGC_PACA_CA_metaMatrix_duct_surv.Rdata")
library(AnnoProbe)
probe2gene=unique(ICGC_PACA_CA_metaMatrix_duct_T[,c("icgc_sample_id","icgc_donor_id.x")])#soft
colnames(probe2gene)=c("probe_id","symbol")

sampletodonor <- function(probes_expr){
  probes_expr=as.data.frame(t(probes_expr))
  genes_expr <- filterEM(probes_expr,probe2gene)
  genes_expr=as.data.frame(t(genes_expr))
  return(genes_expr)
}
icgc_CA_tpm_immune_dornor=sampletodonor(icgc_CA_tpm_immune)

icgc_CA_tpm_immune_dornor_lose3=icgc_CA_tpm_immune_dornor[-c(66,67,68),]
save(icgc_CA_tpm_immune_dornor,icgc_CA_tpm_immune_dornor_lose3,file="ppp_icgc_ca_dornor_lose3.Rdata")

#95
# (ICGC_PACA_CA_seq_log_filter_dornor,ICGC_PACA_CA_seq_norm_dornor,ICGC_PACA_CA_seq_raw_dornor,ICGC_PACA_CA_metaMatrix_duct_surv,ICGC_PACA_CA_metaMatrix_duct_T)
