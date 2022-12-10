rm(list=ls())
options(stringsAsFactors = F)
gc()
library(AnnoProbe)
library(oligo)#不要同时和library(affy)一起load
library(pd.hg.u133.plus.2)
library(readxl)
# library(limma)
# library(hpgltools)#detach("package:limma")#不能一起用
# library(lumi)
library(GEOquery)
# kkk=exprs(affyRaw)
# kkk=normalizeQuantiles(kkk)
# kkk=log2(kkk+1)
# boxplot(kkk,las=2)
#gset=AnnoProbe::geoChina('GSE62254')
load("./GSE62254_eSet.Rdata")#log不平，但像是弄过了
# check the ExpressionSet
eSet=gset[[1]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)#原矩阵有错
head(probes_expr[,1:4])
# probes_expr=log2(probes_expr+1)
# boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])

# if(eSet@annotation=="GPL570"){
system("tar -xvf ./GSE62254/GSE62254_RAW.tar -C ./GSE62254/")
celFiles <- list.celfiles('./GSE62254/',full.name=TRUE,listGzipped = T)
affyRaw <- read.celfiles(celFiles)
# 提取矩阵并做normalization 
#rma 自带 背景校正方法：rma，标准化方法：quantile，汇总方法：medianpolish
#不需要重复quantile，当然limma要干嘛随便
eset <- rma(affyRaw)#不存在直接可以复现的quantile，因为这是用C语言写的
# 检查数据
probes_expr_rma=exprs(eset)
dim(probes_expr_rma)
#
phenoDat_rma <- pData(eset)
head(phenoDat_rma[,1:4])
#
#read_xls处理预后信息
GSE62254_outcome=read.csv("./GSE62254/GSE62254_outcome.csv",header = T)
GSE62254_outcome=as.data.frame(GSE62254_outcome)
#
pdf("GSE62254_boxplot.pdf",width=100)
par(mfrow = c(2,1));
boxplot(probes_expr,las=2)
boxplot(probes_expr_rma,las=2)
dev.off()
#选择 注意每次要改 开始修剪了!!!!!!!!!
probes_expr=probes_expr_rma
colnames(probes_expr)=substr(colnames(probes_expr),1,10)
phenoDat=merge(phenoDat,GSE62254_outcome,by.x="patient:ch1",by.y="Tumor.ID")
rownames(phenoDat)=phenoDat$geo_accession
## check GPL and annotate the probes to genes.
(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl=gpl,type = "bioc")#soft的合并性差落后，建议首选bioc
head(probe2gene)
genes_expr <- filterEM(probes_expr,probe2gene)
head(genes_expr)
# 修剪+median !!!!!!!!!!!!!!!!!!!!!!!!
dim(phenoDat)
phenoDat$time=as.numeric(phenoDat$time)
phenoDat=phenoDat[phenoDat$time>=1,]
dim(phenoDat)
genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
dim(genes_expr_with_clinical)
#
library(DescTools)
genes_expr_mad_GSE62254=t(RobScale(t(genes_expr_with_clinical)))
genes_expr_mean_GSE62254=t(scale(t(genes_expr_with_clinical)))
probes_expr_GSE62254=probes_expr
phenoDat_GSE62254=phenoDat
save(genes_expr_mad_GSE62254,genes_expr_mean_GSE62254,probes_expr_GSE62254,phenoDat_GSE62254,file="GSE62254_after_bioc.Rdata")
write.csv(phenoDat_GSE62254,file="phenoDat_GSE62254.csv",quote = T)
pdf("GSE62254_mad.pdf",width = 100)
par(mfrow = c(2,1));
boxplot(genes_expr_mad)
genes_expr_mean=t(scale(t(genes_expr)))
boxplot(genes_expr_mean)
dev.off()
