rm(list=ls())
options(stringsAsFactors = F)
gc()
library(oligo)#不要同时和library(affy)一起load
library(pd.hg.u133.plus.2)
library(readxl)
# library(limma)
library(hpgltools)#detach("package:limma")#不能一起用
# library(lumi)

# kkk=exprs(affyRaw)
# kkk=normalizeQuantiles(kkk)
# kkk=log2(kkk+1)
# boxplot(kkk,las=2)
#gset=AnnoProbe::geoChina('GSE62254')
load("./GSE62254_eSet.Rdata")#log不平，但像是弄过了
# check the ExpressionSet
eSet=gset[[1]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)
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
eset <- rma(affyRaw,normalize=F)#不存在直接可以复现的quantile，因为这是用C语言写的
# 检查数据
probes_expr_rma=exprs(eset)
dim(probes_expr_rma)
#
phenoDat_rma <- pData(eset)
head(phenoDat_rma[,1:4])
#
GSE62254_outcome=read_xls("./GSE62254/GSE62254_outcome.xls",col_names = TRUE)
GSE62254_outcome=as.data.frame(GSE62254_outcome)
GSE62254_outcome=GSE62254_outcome[-which(substr(GSE62254_outcome$`GSM ID`,1,1)=="*"),]
pdf("GSE62254_boxplot.pdf",width=100)
par(mfrow = c(2,1));
boxplot(probes_expr,las=2)
boxplot(probes_expr_rma,las=2)
dev.off()
save("GSE62254_after.Rdata",probes_expr,phenoDat,probes_expr_rma,_rma,GSE62254_outcome)
# }

normalize.method()
bgcorrect.methods()
showMethods("normalize")
normalize
selectMethod("normalize", "AffyBatch")

p1=normalize(probes_expr, method = )
p2=normalize(probes_expr_rma)#*

