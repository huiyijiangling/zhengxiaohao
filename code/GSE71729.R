#还没改呢 20210309这个地方要改的很多啊。。。
#quantile normalize
rm(list=ls())
options(stringsAsFactors = F)
gc()
library(AnnoProbe)
library(oligo)#不要同时和library(affy)一起load
library(pd.hg.u133.plus.2)
library(readxl)
library(DescTools)
library(miRBaseConverter)
library(tidyr)
library(stringr)
# library(limma)
# library(hpgltools)#detach("package:limma")#不能一起用
# library(lumi)
library(GEOquery)
# kkk=exprs(affyRaw)
# kkk=normalizeQuantiles(kkk)
# kkk=log2(kkk+1)
# boxplot(kkk,las=2)
#gset=AnnoProbe::geoChina('GSE71729')
load("./GSE71729_eSet.Rdata")#log不平，但像是弄过了
# check the ExpressionSet
eSet=gset[[1]]#eSet=gset[[2]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)#原矩阵,途经过滤1次
head(probes_expr[,1:4])
# probes_expr=log2(probes_expr+1)
# boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
if(eSet@annotation=="GPL570"){
  system("tar -xvf ./GSE71729/GSE71729_RAW.tar -C ./GSE71729/")
  celFiles <- list.celfiles('./GSE71729/',full.name=TRUE,listGzipped = T)
  affyRaw <- read.celfiles(celFiles)
  # 提取矩阵并做normalization 
  #rma 自带 背景校正方法：rma，标准化方法：quantile，汇总方法：medianpolish
  #不需要重复quantile，当然limma要干嘛随便
  eset <- rma(affyRaw)#不存在直接可以复现的quantile，因为这是用C语言写的
  # 检查数据
  probes_expr_rma=exprs(eset)
  dim(probes_expr_rma)
  #
  pdf("GSE71729_boxplot.pdf",width=100)
  par(mfrow = c(2,1));
  boxplot(probes_expr,las=2)
  boxplot(probes_expr_rma,las=2)
  dev.off()
  #
  phenoDat_rma <- pData(eset)
  head(phenoDat_rma[,1:4])
  # 5 P/A过滤
  # 目的是去除“不表达”的基因/探针数据，使用paCalls函数，选取p值小于0.05的探针：
  xpa <- paCalls(affyRaw)
  View(xpa)#要看返回的是一列还是两列
  AP <- apply(xpa$p, 1, function(x) any(x < 0.05))
  xids <- names(AP[AP])#xids <- as.numeric(names(AP[AP]))
  head(xids)
  ## [1] 11 14 17 20 28 42
  
  # 计算后发现P/A结果和表达量分析结果所用的探针名称是不一样的，需要使用探针信息进行转换。getProbeInfo函数可以获取相关信息：
  pinfo <- getProbeInfo(affyRaw)
  head(pinfo)
  # fid   man_fsetid
  # 1 1170    236072_at
  # 2 1171    206490_at
  # 3 1172 1568408_x_at
  # 4 1173    235895_at
  # 5 1174  211002_s_at
  # 6 1175    239356_at
  # 
  # 两类id转换后进行筛选：
  
  fids <- pinfo[pinfo$man_fsetid %in% xids, 2] #fids <- pinfo[pinfo$fid %in% xids, 2]
  head(fids)
  ## [1] "13402401" "13363588" "13492405" "13488138" "13349578" "13529865"
  nrow(probes_expr_rma)
  ## [1] 38408
  probes_expr_rma <- probes_expr_rma[rownames(probes_expr_rma) %in% fids, ]
  nrow(probes_expr_rma)
  ## [1] 35697
  # 第一个语句筛选出表达基因的探针ID，第四个语句用筛选出的ID过滤掉不表达的基因/探针。
  #选择 注意每次要改 开始修剪了!!!!!!!!!
  probes_expr=probes_expr_rma
}
# colnames(probes_expr)=substr(colnames(probes_expr),1,10)
#
#read_xls处理预后信息 没有
# phenoDat_outcome=read.csv("./GSE71729/GSE71729_outcome.csv",header = T)
# phenoDat_outcome=as.data.frame(phenoDat_outcome)
#
# phenoDat=merge(phenoDat,phenoDat_outcome,by.x="patient:ch1",by.y="Tumor.ID")
# rownames(phenoDat)=phenoDat$geo_accession
## check GPL and annotate the probes to genes.
if(F){
(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl=gpl,type = "soft")#soft的合并性差落后，建议首选bioc,但如果结合了更新就不一样了,例如ta
head(probe2gene)
}
#为了图方便我用了直接下载的矩阵!!!!!
probe2gene=as.data.frame(probes_expr)
probe2gene$symbol=rownames(probes_expr)
probe2gene$ID=rownames(probes_expr)
#可能不行啊（symbol滤过了少了很多啊）
# source("./updateName.R")
# forthetree=updateName(probe2gene$symbol)
# table(forthetree$gene_biotype)
# probe2gene=merge(forthetree,probe2gene,by.x="ALIAS",by.y="symbol")
# probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
# probe2gene=unique(probe2gene[,c("probe_id","ENSEMBL")])
# colnames(probe2gene)=c("probe_id","symbol")
# genes_expr <- filterEM(probes_expr,probe2gene)
# head(genes_expr)
source("./updateName.R")# ncbi与biomart相比多map了200，等于是，ncbi》biomart》org.hs
forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol)
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
# probe2gene=unique(probe2gene[,c("probe_id","ENSEMBL")])#bioc
probe2gene=unique(probe2gene[,c("ID","ENSEMBL")])#soft
# probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
colnames(probe2gene)=c("probe_id","symbol")
genes_expr <- filterEM(probes_expr,probe2gene)
# 修剪+median !!!!!!!!!!!!!!!!!!!!!!!!
dim(phenoDat)
# phenoDat$time=as.numeric(phenoDat$time)
# phenoDat=phenoDat[phenoDat$time>=1,]
# dim(phenoDat)
phenoDat=phenoDat[phenoDat$characteristics_ch2.1=="tissue type: Primary",]
phenoDat=phenoDat[!is.na(phenoDat$`death_event_1death_0censor:ch2`),]

phenoDat$OS=as.numeric(phenoDat$`death_event_1death_0censor:ch2`)
phenoDat$OS.time=as.numeric(phenoDat$`survival_months:ch2`)
phenoDat_GSE71729_rna_125_keepall=phenoDat

genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
genes_expr_mean_GSE71729_rna_125=t(scale(t(genes_expr_with_clinical)))

phenoDat=subset(phenoDat,OS.time>=1)

genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
dim(genes_expr_with_clinical)
genes_expr_GSE71729_rna=genes_expr_with_clinical
phenoDat_GSE71729_rna=phenoDat
probes_expr_GSE71729_rna=probes_expr
genes_expr_mad_GSE71729_rna=t(RobScale(t(genes_expr_with_clinical)))
genes_expr_mean_GSE71729_rna=t(scale(t(genes_expr_with_clinical)))


#我觉得不用删掉，只是缺东西而已。全删除就是97个、


#
save(phenoDat_GSE71729_rna_125_keepall,genes_expr_mean_GSE71729_rna_125,phenoDat_GSE71729_rna,probes_expr_GSE71729_rna,genes_expr_mad_GSE71729_rna,genes_expr_mean_GSE71729_rna,genes_expr_GSE71729_rna,
     file="GSE71729_after_soft.Rdata")
#(file="GSE71729_after_bioc.Rdata"))
write.csv(phenoDat_GSE71729_rna_125_keepall,file="phenoDat_GSE71729_rna_125_keepall.csv",quote = T)
write.csv(phenoDat_GSE71729_rna,file="phenoDat_GSE71729_rna.csv",quote = T)
pdf("GSE71729_mad.pdf",width = 100)
par(mfrow = c(2,1));
boxplot(probes_expr)
genes_expr_mean=t(scale(t(genes_expr)))
boxplot(genes_expr_mean)
dev.off()

