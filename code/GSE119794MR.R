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
#gset=AnnoProbe::geoChina('GSE119794')
load("./GSE119794_eSet.Rdata")#log不平，但像是弄过了
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
if(T){
  # gset=getGEO('GSE119794',getGPL = F)
  # getGEOSuppFiles('GSE119794')
  
  #windows powershell
  # 更改文件权限
  system("tar -xvf ./GSE119794/GSE119794_RAW.tar -C ./GSE119794/")
  
  #或者        filenames <- file.path(path, metadata$file_id, metadata$file_name, fsep = .Platform$file.sep)
  filenames_M <- list.files('./GSE119794/',pattern = "M.txt.gz",full.name=TRUE)
  mirCount <- do.call("cbind", lapply(filenames_M, function(fl) 
    read.delim(fl)[,2]))
  rownames(mirCount) <- read.delim(filenames_M[1])[,1]
  fn=unlist(lapply(strsplit(unlist(lapply(strsplit(filenames_M,'/'), function(x)x[3])),"_"),function(x)x[1]))
  colnames(mirCount)=phenoDat[match(fn,phenoDat$geo_accession),2]
  
  filenames_R <- list.files('./GSE119794/',pattern = "R.txt.gz",full.name=TRUE)
  rnaCount <- do.call("cbind", lapply(filenames_R, function(fl) 
    read.delim(fl)[,2]))
  rownames(rnaCount) <- read.delim(filenames_R[1])[,1]
  fn=unlist(lapply(strsplit(unlist(lapply(strsplit(filenames_R,'/'), function(x)x[3])),"_"),function(x)x[1]))
  colnames(rnaCount)=phenoDat[match(fn,phenoDat$geo_accession),2]
  probes_expr_GSE119794_mir=mirCount
  probes_expr_GSE119794_rna=rnaCount
  phenoDat_GSE119794=phenoDat
  save(probes_expr_GSE119794_mir,probes_expr_GSE119794_rna,phenoDat_GSE119794,file ="GSE119794_before.Rdata" )
  write.csv(phenoDat_GSE119794,file="phenoDat_GSE119794.csv",quote = T)
  #https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119794/matrix/GSE119794_series_matrix.txt.gz
  #https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119794/matrix
  #https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119794/suppl/ 如果没上传就放在这
  #有时候原来的矩阵是没有上传的 那就可以去getGEO   
  #getGEOfile('GSE119794') 这个没有用
}
if(eSet@annotation=="GPL570"){
  system("tar -xvf ./GSE119794/GSE119794_RAW.tar -C ./GSE119794/")
  celFiles <- list.celfiles('./GSE119794/',full.name=TRUE,listGzipped = T)
  affyRaw <- read.celfiles(celFiles)
  # 提取矩阵并做normalization 
  #rma 自带 背景校正方法：rma，标准化方法：quantile，汇总方法：medianpolish
  #不需要重复quantile，当然limma要干嘛随便
  eset <- rma(affyRaw)#不存在直接可以复现的quantile，因为这是用C语言写的
  # 检查数据
  probes_expr_rma=exprs(eset)
  dim(probes_expr_rma)
  #
  pdf("GSE119794_boxplot.pdf",width=100)
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
#read_xls处理预后信息
# phenoDat_outcome=read.csv("./GSE119794/GSE119794_outcome.csv",header = T)
# phenoDat_outcome=as.data.frame(phenoDat_outcome)
#
# phenoDat=merge(phenoDat,phenoDat_outcome,by.x="patient:ch1",by.y="Tumor.ID")
# rownames(phenoDat)=phenoDat$geo_accession
## check GPL and annotate the probes to genes.
(gpl=eSet@annotation)
checkGPL(gpl)
#illumina2000在里面没有额，只好自己弄了
load("GSE119794_before.Rdata")
probe2gene=data.frame(probe_id=rownames(probes_expr_GSE119794_rna),symbol=rownames(probes_expr_GSE119794_rna))
source("./updateName.R")
forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol,gene.info.file = "Homo_sapiens.gene_info.gz")
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
if(F){
forthetree=updateName(probe2gene$symbol)
table(forthetree$gene_biotype)
probe2gene=merge(forthetree,probe2gene,by.x="ALIAS",by.y="symbol",all.y=T)
}
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
probe2gene=unique(probe2gene[,c("probe_id","ENSEMBL")])
# probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
colnames(probe2gene)=c("probe_id","symbol")
probes_expr=probes_expr_GSE119794_rna
genes_expr <- filterEM(probes_expr,probe2gene)
genes_expr_GSE119794_rna=genes_expr
# probes_expr_GSE119794_rna["MIR4435-1HG",]
# genes_expr["ENSG00000172965",]


(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
GPL=getGEO(gpl,AnnotGPL = F,GSEMatrix = F, getGPL =T,parseCharacteristics =T)
GPL=GPL@dataTable@table
#只好自己造了
GPL=data.frame(rownames(probes_expr_GSE119794_mir),rownames(probes_expr_GSE119794_mir))
colnames(GPL)=c("ID","miRNA_ID")
miRNANames=rownames(probes_expr_GSE119794_mir)
version=checkMiRNAVersion(miRNANames, verbose = FALSE)
version#12
tbb= miRNA_NameToAccession(miRNANames,version = version)#只有一个了
tbb=unique(tbb)
tbb2=miRNA_AccessionToName(tbb$Accession,targetVersion = "v21")
tbb3=merge(tbb,tbb2,by="Accession",all=T)
probe2gene=merge(GPL,tbb3,by.x="miRNA_ID",by.y=paste0("miRNAName_",version,sep=""),all=T)
probe2gene=probe2gene[!is.na(probe2gene$TargetName),]
probe2gene=probe2gene[(str_split(probe2gene$TargetName,'-',simplify = T)[,1]) == "hsa",]#不筛掉人类
# GPL=GPL[GPL$CONTROL_TYPE=="FALSE",]
probe2gene=unique(probe2gene[,c("ID","TargetName")])
colnames(probe2gene)=c("probe_id","symbol")
table(duplicated(probe2gene$probe_id))
table(duplicated(probe2gene$symbol))
probes_expr=probes_expr_GSE119794_mir
# printGPLInfo(gpl)
# probe2gene=idmap(gpl=gpl,type = "bioc")#soft的合并性差落后，建议首选bioc
# head(probe2gene)
genes_expr <- filterEM(probes_expr,probe2gene)
genes_expr_GSE119794_mir=genes_expr
# head(genes_expr)
# 修剪+median !!!!!!!!!!!!!!!!!!!!!!!!
# dim(phenoDat)
# phenoDat$time=as.numeric(phenoDat$time)
# phenoDat=phenoDat[phenoDat$time>=1,]
# dim(phenoDat)
# genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
# dim(genes_expr_with_clinical)
save(probes_expr_GSE119794_mir,genes_expr_GSE119794_mir,
     probes_expr_GSE119794_rna,genes_expr_GSE119794_rna,
     phenoDat_GSE119794,file ="GSE119794_before.Rdata")


###########################
if(F){
library(DescTools)
#
genes_expr_mad_GSE119794_rna=t(RobScale(t(rnaCount)))
genes_expr_mean_GSE119794_rna=t(scale(t(rnaCount)))


genes_expr_mad_GSE119794_mir=t(RobScale(t(mirCount)))
genes_expr_mean_GSE119794_mir=t(scale(t(mirCount)))



save(genes_expr_mad_GSE119794,genes_expr_mean_GSE119794,probes_expr_GSE119794,phenoDat_GSE119794,file="GSE119794_after_bioc.Rdata")
write.csv(phenoDat_GSE119794,file="phenoDat_GSE119794.csv",quote = T)
pdf("GSE119794_mad.pdf",width = 100)
par(mfrow = c(2,1));
boxplot(genes_expr_mad)
genes_expr_mean=t(scale(t(genes_expr)))
boxplot(genes_expr_mean)
dev.off()
}