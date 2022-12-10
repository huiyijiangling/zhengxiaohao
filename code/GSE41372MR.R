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
#gset=AnnoProbe::geoChina('GSE41372')
load("./GSE41372_eSet.Rdata")#log不平，但像是弄过了
# check the ExpressionSet
eSet=gset[[2]]#eSet=gset[[2]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)#原矩阵有错
head(probes_expr[,1:4])
# probes_expr=log2(probes_expr+1)
# boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
if(eSet@annotation=="GPL570"){
  system("tar -xvf ./GSE41372/GSE41372_RAW.tar -C ./GSE41372/")
  celFiles <- list.celfiles('./GSE41372/',full.name=TRUE,listGzipped = T)
  affyRaw <- read.celfiles(celFiles)
  # 提取矩阵并做normalization 
  #rma 自带 背景校正方法：rma，标准化方法：quantile，汇总方法：medianpolish
  #不需要重复quantile，当然limma要干嘛随便
  eset <- rma(affyRaw)#不存在直接可以复现的quantile，因为这是用C语言写的
  # 检查数据
  probes_expr_rma=exprs(eset)
  dim(probes_expr_rma)
  #
  pdf("GSE41372_boxplot.pdf",width=100)
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
# phenoDat_outcome=read.csv("./GSE41372/GSE41372_outcome.csv",header = T)
# phenoDat_outcome=as.data.frame(phenoDat_outcome)
#
# phenoDat=merge(phenoDat,phenoDat_outcome,by.x="patient:ch1",by.y="Tumor.ID")
# rownames(phenoDat)=phenoDat$geo_accession
## check GPL and annotate the probes to genes.
(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl=gpl,type = "soft")#soft的合并性差落后，建议首选bioc,但如果结合了更新就不一样了,例如ta
head(probe2gene)
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
# phenoDat=phenoDat[phenoDat$source_name_ch1 %in% c("ductal adenocarcinoma","non-neoplastic pancreas"),]
genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
dim(genes_expr_with_clinical)
genes_expr_GSE41372_rna=genes_expr_with_clinical
phenoDat_GSE41372_rna=phenoDat
probes_expr_GSE41372_rna=probes_expr
genes_expr_mad_GSE41372_rna=t(RobScale(t(genes_expr_with_clinical)))
genes_expr_mean_GSE41372_rna=t(scale(t(genes_expr_with_clinical)))
#
eSet=gset[[1]]#eSet=gset[[2]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)#原矩阵有错
head(probes_expr[,1:4])
# probes_expr=log2(probes_expr+1)
# boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])

(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
GPL=getGEO(gpl,AnnotGPL = F,GSEMatrix = F, getGPL =T,parseCharacteristics =T)
GPL=GPL@dataTable@table
# GPL=GPL[GPL$CONTROL_TYPE=="FALSE",]
# probe2gene=GPL[,c("ID","GENE_SYMBOL")]
# colnames(probe2gene)=c("probe_id","symbol")
GPL[GPL==""] <- NA
table(GPL[!is.na(GPL$miRNA_LIST) & is.na(GPL$Accession),])#要有一起有
miRNANames=GPL[is.na(GPL$miRNA_LIST),]$miRNA_LIST
# miRNANames=GPL[is.na(GPL$Accession),]$SPOT_ID
taa = miRNAVersionConvert(miRNANames,targetVersion = "v21",exact = TRUE)#这是用来找NA的
table(taa$TargetName)
# 分列名称和accession号
# 注意一定不要把ID拿去裂项，就算错了一起错还是能够配上
GPL=separate_rows(GPL,Accession,sep = "[^-*[:alnum:].]+")
GPL=separate_rows(GPL,miRNA_LIST,sep = "[^-*[:alnum:].]+")
GPL=unique(GPL)
miRNANames=GPL[!is.na(GPL$miRNA_LIST),]$miRNA_LIST
version=checkMiRNAVersion(miRNANames, verbose = FALSE)
version#12
tbb= miRNA_NameToAccession(miRNANames,version = version)#只有一个了
tbb=unique(tbb)
tbb2=miRNA_AccessionToName(tbb$Accession,targetVersion = "v21")
tbb3=merge(tbb,tbb2,by="Accession",all=T)
probe2gene=merge(GPL,tbb3,by.x="miRNA_LIST",by.y=paste0("miRNAName_",version,sep=""),all=T)
probe2gene=probe2gene[!is.na(probe2gene$TargetName),]
probe2gene=probe2gene[(str_split(probe2gene$TargetName,'-',simplify = T)[,1]) == "hsa",]#不筛掉人类
# GPL=GPL[GPL$CONTROL_TYPE=="FALSE",]
probe2gene=unique(probe2gene[,c("ID","TargetName")])
colnames(probe2gene)=c("probe_id","symbol")
table(duplicated(probe2gene$probe_id))
table(duplicated(probe2gene$symbol))
probe2gene=unlist(lapply(split(probe2gene$symbol,f=probe2gene$probe_id) ,function(x) paste0(x,collapse="/")))
probe2gene <- as.data.frame(probe2gene)
colnames(probe2gene)="symbol"
probe2gene$probe_id=rownames(probe2gene)
orderorder = c("probe_id","symbol")
probe2gene = probe2gene[orderorder]
# probe2gene=idmap(gpl=gpl,type = "bioc")#soft的合并性差落后，建议首选bioc
# head(probe2gene)
genes_expr <- filterEM(probes_expr,probe2gene)
# head(genes_expr)
# 修剪+median !!!!!!!!!!!!!!!!!!!!!!!!
dim(phenoDat)
# phenoDat$time=as.numeric(phenoDat$time)
# phenoDat=phenoDat[phenoDat$time>=1,]
# dim(phenoDat)
# phenoDat=phenoDat[phenoDat$source_name_ch1 %in% c("ductal adenocarcinoma","non-neoplastic pancreas"),]
genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
dim(genes_expr_with_clinical)
phenoDat_GSE41372_mir=phenoDat
probes_expr_GSE41372_mir=probes_expr
genes_expr_GSE41372_mir=genes_expr_with_clinical
genes_expr_mad_GSE41372_mir=t(RobScale(t(genes_expr_with_clinical)))#可以进行NA检验
genes_expr_mean_GSE41372_mir=t(scale(t(genes_expr_with_clinical)))

#手动检验错误
#手动检验错误
#手动检验错误
#手动检验错误
#手动检验错误
#手动检验错误
phenoDat_GSE41372_rna["GSM1015737","title"]="pancreatic ductal adenocarcinoma - Patient Sample 09"

# 合并为了留下统一用户
phenoDat_GSE41372_mir=phenoDat_GSE41372_mir[phenoDat_GSE41372_mir$title %in% phenoDat_GSE41372_rna$title,]
genes_expr_GSE41372_mir=genes_expr_GSE41372_mir[,which(colnames(genes_expr_GSE41372_mir) %in% phenoDat_GSE41372_mir$geo_accession)]
phenoDat_GSE41372_rna=phenoDat_GSE41372_rna[phenoDat_GSE41372_rna$title %in% phenoDat_GSE41372_mir$title,]
genes_expr_GSE41372_rna=genes_expr_GSE41372_rna[,which(colnames(genes_expr_GSE41372_rna) %in% phenoDat_GSE41372_rna$geo_accession)]

#
save(phenoDat_GSE41372_mir,probes_expr_GSE41372_mir,genes_expr_mad_GSE41372_mir,genes_expr_mean_GSE41372_mir,genes_expr_GSE41372_mir,
     phenoDat_GSE41372_rna,probes_expr_GSE41372_rna,genes_expr_mad_GSE41372_rna,genes_expr_mean_GSE41372_rna,genes_expr_GSE41372_rna,
     file="GSE41372_after_soft.Rdata")
#(file="GSE41372_after_bioc.Rdata"))
# write.csv(phenoDat_GSE41372,file="phenoDat_GSE41372.csv",quote = T)
pdf("GSE41372_mad.pdf",width = 100)
par(mfrow = c(2,1));
boxplot(genes_expr_mad)
genes_expr_mean=t(scale(t(genes_expr)))
boxplot(genes_expr_mean)
dev.off()