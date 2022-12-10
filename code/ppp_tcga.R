###################################################################################
rm(list=ls())
options(stringsAsFactors = F)
gc()
#######################################1 数据准备
load("TcgaTargetGTEX_phenotype_Pancreas.Rdata")
load("rnaCountslog21_Pancreas.Rdata")
load("rnatpmlog2001_Pancreas.Rdata")
rownames(TcgaTargetGTEX_phenotype_Pancreas[TcgaTargetGTEX_phenotype_Pancreas$sample_type=="WRONG",])#"TCGA-HZ-A9TJ-06"
TcgaTargetGTEX_phenotype_Pancreas624=TcgaTargetGTEX_phenotype_Pancreas[-which(rownames(TcgaTargetGTEX_phenotype_Pancreas) %in% c("TCGA-HZ-A9TJ-06")),]
rownames(rnatpmlog2001_Pancreas)=substr(rownames(rnatpmlog2001_Pancreas),1,15)
rownames(rnaCountslog21_Pancreas)=substr(rownames(rnaCountslog21_Pancreas),1,15)
rnatpmlog2001_Pancreas=subset(rnatpmlog2001_Pancreas,select=rownames(TcgaTargetGTEX_phenotype_Pancreas624))
rnaCountslog21_Pancreas=subset(rnaCountslog21_Pancreas,select=rownames(TcgaTargetGTEX_phenotype_Pancreas624))
rnatpm_Pancreas=2^rnatpmlog2001_Pancreas-0.001#625
rnaCounts_Pancreas=2^rnaCountslog21_Pancreas-1#624
# write.csv(rnaCounts_Pancreas,file="paad.csv",quote=T)
# write.csv(rnaCounts,file="60843.csv",quote=T)

# ######################################2 基因长度计算
# library("GenomicFeatures")
# ## 导入gff3文件
# txdb <- makeTxDbFromGFF("gencode.v36.annotation.gtf.gz",format="gtf")
# ## 获取外显子位置
# exons_gene <- exonsBy(txdb, by = "gene")
# ## 去除外显子重叠部分，计算外显子长度
# exons_gene_len <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
# exons_gene_len2 <- do.call(rbind,lapply(exons_gene_len, data.frame))
# exons_gene_len2$gene_id=stringr::str_split(rownames(exons_gene_len2),"\\.",simplify=T)[,1]
# exons_gene_len3=exons_gene_len2[!duplicated(exons_gene_len2$gene_id),]
# rownames(exons_gene_len3)=exons_gene_len3$gene_id
# exons_gene_len3$eff_length=exons_gene_len3$X..i..
# exons_gene_len3=exons_gene_len3[,-1]
# exons_gene_len=exons_gene_len3
# save(exons_gene_len,file="exonlength.Rdata")
# ###################################### 3count to tpm
# load("exonlengthv22.Rdata")
# exons_gene_len[1:3,]#还正常，只不过版本确实有点不一样
# # ENSG00000000003 ENSG00000000003       4536
# # ENSG00000000005 ENSG00000000005       1476
# # ENSG00000000419 ENSG00000000419       1207
# # gene_id eff_length v22 也不行啊
# # ENSG00000000003 ENSG00000000003       4535
# # ENSG00000000005 ENSG00000000005       1610
# # ENSG00000000419 ENSG00000000419       1207
# gen <- intersect(rownames(rnaCounts_Pancreas),rownames(exons_gene_len))
# rnaCounts_Pancreas_tpm=rnaCounts_Pancreas[gen,]
# eff_length=exons_gene_len[gen,]
# countToTpm <- function(counts, effLen)
# {
#   rate <- log(counts) - log(effLen)
#   denom <- log(sum(exp(rate)))
#   exp(rate - denom + log(1e6))
# }
# countToFpkm <- function(counts, effLen)
# {
#   N <- sum(counts)
#   exp( log(counts) + log(1e9) - log(effLen) - log(N) )
# }
# 
# rnaCounts_Pancreas_tpm <- apply(rnaCounts_Pancreas_tpm,2,function(x) countToTpm(x,eff_length$eff_length))
# colSums(rnaCounts_Pancreas_tpm)
# ###################################### 4check
# median(rnaCounts_Pancreas_tpm["ENSG00000143947",])
# median(rnatpm_Pancreas["ENSG00000143947",])
# log2(rnaCounts_Pancreas_tpm["ENSG00000143947",]+0.001)
# rnatpmlog2001_Pancreas["ENSG00000143947",]


# dataset=rnatpm_Pancreas
# dataset[is.na(dataset)] <- 0
# load("C:/Users/zxh/Desktop/R/gtf/exonlength.Rdata")


# gen <- intersect(rownames(dataset),rownames(exons_gene_len))
# dataset_tpm=dataset[gen,]
# eff_length=exons_gene_len[gen,]
# countToTpm <- function(counts, effLen)
# {
#   rate <- log(counts) - log(effLen)
#   denom <- log(sum(exp(rate)))
#   exp(rate - denom + log(1e6))
# }
# countToFpkm <- function(counts, effLen)
# {
#   N <- sum(counts)
#   exp( log(counts) + log(1e9) - log(effLen) - log(N) )
# }
# 
# dataset_tpm <- apply(dataset_tpm,2,function(x) countToTpm(x,eff_length$eff_length))
# colSums(dataset_tpm)
###################################### 4check
median(dataset_tpm["ENSG00000143947",])
median(rnatpm_Pancreas["ENSG00000143947",])
log2(dataset_tpm["ENSG00000143947",]+0.001)
rnatpmlog2001_Pancreas["ENSG00000143947",]

dataset=dataset_tpm


# library(AnnoProbe)
# probe2gene=data.frame(probe_id=rownames(dataset),symbol=rownames(dataset))
# source("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/updateName.R")
# forthetree=ensg2SymbolUsingNCBImodify(alias=probe2gene$symbol,gene.info.file = "Homo_sapiens.gene_info.gz")
# probe2gene=cbind(probe2gene,forthetree)
# table(is.na(probe2gene$ENSEMBL))
# if(F){
#   forthetree=updateName(probe2gene$symbol)
#   table(forthetree$gene_biotype)
#   probe2gene=merge(forthetree,probe2gene,by.x="ALIAS",by.y="symbol",all.y=T)
# }
# probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
# probe2gene=unique(probe2gene[,c("probe_id","Symbol")])
# # probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
# colnames(probe2gene)=c("probe_id","symbol")
# probes_expr=dataset
# genes_expr <- filterEM(probes_expr,probe2gene)
# ICGC_PACA_CA_seq_raw_symbol_tpm=genes_expr
# write.csv(ICGC_PACA_CA_seq_raw_symbol_tpm,"ICGC_PACA_CA_seq_raw_symbol_tpm.csv",quote = T)


#载入R包；
# library(ggcorrplot)
library(ggplot2)
#直接快速绘制整个相关性热图；

# quickcor(mtcars, cluster = TRUE,cor.test = TRUE) +
#   geom_colour() +
#   geom_mark(size=3,color="white",fontface=1)+
#   scale_fill_gradientn(colours = c("#77C034","white" ,"#C388FE"))+
#   geom_panel_grid(colour = "white",size = 1)

TCGA_immune=readxl::read_xlsx("C:/Users/zxh/Desktop/R/ftoimmune/xCell_TCGA_RSEM.xlsx")
TCGA_immune=as.data.frame(TCGA_immune)
rownames(TCGA_immune)=TCGA_immune[,1]
TCGA_immune=TCGA_immune[,-1]
TcgaTargetGTEX_phenotype_Pancreas624T=subset(TcgaTargetGTEX_phenotype_Pancreas624,rownames(TcgaTargetGTEX_phenotype_Pancreas624)%in%colnames(TCGA_immune))
TCGA_immune=subset(TCGA_immune,select=rownames(TcgaTargetGTEX_phenotype_Pancreas624T))
cor=as.data.frame(t(TCGA_immune))
rnatpm_Pancreast=as.data.frame(t(rnatpm_Pancreas))
cor$sample=rownames(cor)
rnatpm_Pancreast$sample=rownames(rnatpm_Pancreast)
cor=merge(cor,rnatpm_Pancreast,by="sample")
cor
cor=subset(cor,select=c("ENSG00000140718",rownames(TCGA_immune)))
library(tidyverse)
library(corrplot)
# p.mat=res1$p,

# cor=t(cor)
cor=as.data.frame(cor)
colnames(cor)=c("FTO",rownames(TCGA_immune))
# cor=as.numercoric(cor)
res1 <- cor.mtest(cor, conf.level = .95, method = "spearman")
cor(cor) %>%   corrplot(method = "square",
                        type = "lower",tl.srt = 45,tl.col = "black")

pdf("CoR.pdf",width = 10,height=10)
pp1<- cor(cor, method ="spearman") %>%  corrplot(type="lower",method = "square",
                                                 insig="label_sig",p.mat = res1$p,outline = "white",
                                                 sig.level=c(.05),  ,#.001,.01,
                                                 pch.cex=.9,pch.col="green",tl.srt=45,tl.col="black")
dev.off()
# order = "AOE", cl.pos = "b", tl.pos = "d", tl.srt = 60) 


# pdf("Construction_risk_score.pdf",width = 8,height=8)
# cor(cor) %>%   corrplot(method = "square",type = "lower",tl.srt = 45,tl.col = "black",)
# dev.off()
pp2=res1[["p"]]
pp2=as.data.frame(pp2)
pp1=as.data.frame(pp1)
pp1
ppp_tcga=cbind(pp1[,1],pp2[,1])
rownames(ppp_tcga)=c("FTO",rownames(TCGA_immune))
write.csv(ppp_tcga,"ppp_tcga.csv")
save(ppp_tcga,file="ppp_tcga.Rdata")
###########################################################