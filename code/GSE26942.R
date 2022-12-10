GSE26942=GSE26901+GSE26899
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
#gset=AnnoProbe::geoChina('GSE26942')
load("./GSE26942_eSet.Rdata")#log不平，但像是弄过了
# check the ExpressionSet
eSet=gset[[1]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
probes_expr=log2(probes_expr+1)

## check GPL and annotate the probes to genes.
(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl=gpl,type = "bioc")#soft的合并性差落后，建议首选bioc
head(probe2gene)
genes_expr <- filterEM(probes_expr,probe2gene)
head(genes_expr)
# boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
phenoDat$event=substr(phenoDat$characteristics_ch1.5,7,nchar(phenoDat$characteristics_ch1.5))
phenoDat$time=substr(phenoDat$characteristics_ch1.6,27,nchar(phenoDat$characteristics_ch1.6))
# 修剪+median !!!!!!!!!!!!!!!!!!!!!!!!
dim(phenoDat)
phenoDat$time=as.numeric(phenoDat$time)
phenoDat=phenoDat[phenoDat$time>=1,]
dim(phenoDat)
genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
dim(genes_expr_with_clinical)
#
library(DescTools)
genes_expr_mad_GSE26942=t(RobScale(t(genes_expr_with_clinical)))
genes_expr_mean_GSE26942=t(scale(t(genes_expr_with_clinical)))
probes_expr_GSE26942=probes_expr
phenoDat_GSE26942=phenoDat

save(genes_expr_mad_GSE26942,genes_expr_mean_GSE26942,probes_expr_GSE26942,phenoDat_GSE26942,file="GSE26942_after_bioc.Rdata")

#
if(eSet@annotation=="GPL570"){
  system("tar -xvf ./GSE26942/GSE26942_RAW.tar -C ./GSE26942/")
  celFiles <- list.celfiles('./GSE26942/',full.name=TRUE,listGzipped = T)
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
  GSE26942_outcome=read_xls("./GSE26942/GSE26942_outcome.xls",col_names = TRUE)
  GSE26942_outcome=as.data.frame(GSE26942_outcome)
  GSE26942_outcome=GSE26942_outcome[-which(substr(GSE26942_outcome$`GSM ID`,1,1)=="*"),]
  
  save(file="GSE26942_after.Rdata",probes_expr,phenoDat,probes_expr_rma,phenoDat_rma,GSE26942_outcome)
}
#
pdf("GSE26942_boxplot.pdf",width=100)
par(mfrow = c(2,1));
boxplot(probes_expr,las=2)
boxplot(probes_expr_rma,las=2)
dev.off()
#中间去需要矫正一下combat 76，357,
# mod = model.matrix(~as.factor(sample_type), data=TcgaTargetGTEX_phenotype_stomach624)
# # ,mod = NULL,
# rnaExpr_quant=ComBat(dat=rnaExpr_quant,batch=TcgaTargetGTEX_phenotype_stomach624$`_study`,mod=mod)
# rnaExpr=as.data.frame(rnaExpr_quant)