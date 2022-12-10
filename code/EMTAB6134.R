rm(list=ls())
options(stringsAsFactors = F)
gc()
library(AnnoProbe)
library(oligo)#不要同时和library(affy)一起load
library(pd.hg.u133.plus.2)
library(pd.hg.u219)
library(readxl)
# library(limma)
# library(hpgltools)#detach("package:limma")#不能一起用
# library(lumi)
library(GEOquery)
probes_expr <- read.table(file = "./E-MTAB-6134/ProcessedExpression.tsv",header = T,fileEncoding = "UTF-8")
head(probes_expr[,1:4])
# if(eSet@annotation=="GPL570"){
#   system("tar -xvf ./E-MTAB-6134/E-MTAB-6134_RAW.tar -C ./E-MTAB-6134/")

if(F){
if(T){
  celFiles <- list.celfiles('./E-MTAB-6134/original',full.name=TRUE,listGzipped = T)
  affyRaw <- read.celfiles(celFiles)
  # 提取矩阵并做normalization 
  #rma 自带 背景校正方法：rma，标准化方法：quantile，汇总方法：medianpolish
  #不需要重复quantile，当然limma要干嘛随便
  eset <- rma(affyRaw)#不存在直接可以复现的quantile，因为这是用C语言写的
  # 检查数据
  probes_expr_rma=exprs(eset)
  dim(probes_expr_rma)
  #
  pdf("E-MTAB-6134_boxplot.pdf",width=100)
  par(mfrow = c(2,1));
  boxplot(probes_expr,las=2)
  boxplot(probes_expr_rma,las=2)
  dev.off()
  #
  phenoDat_rma <- pData(eset)
  head(phenoDat_rma[,1:4])
}  
if(F){
  # 5 P/A过滤
  # 目的是去除“不表达”的基因/探针数据，使用paCalls函数，选取p值小于0.05的探针：
  xpa <- oligo::paCalls(affyRaw)
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
  
}

  library(genefilter)
  eset.filter <- nsFilter(eset, require.entrez=F, remove.dupEntrez=F,var.cutoff=0.5) #不使用函数的合并探针功能，50%很好有文献支持
  eset.filter$filter.log #查看每一步筛掉多少探针
  eset.filter <- eset.filter$eset #只留下ExpressionSet对象
  probes_expr_eset.filter=exprs(eset.filter)
  probes_expr_eset.filter=as.data.frame(probes_expr_eset.filter)
  pdf("E-MTAB-6134_probes_expr_eset.filter.pdf",width=100)
  par(mfrow = c(2,1));
  boxplot(probes_expr_eset.filter,las=2)
  dev.off()
  #选择 注意每次要改 开始修剪了!!!!!!!!!
  # probes_expr=probes_expr_rma

  probes_expr=probes_expr[rownames(probes_expr) %in% rownames(probes_expr_eset.filter),]
  
  
  #
  # source("filter_f1000.R")
  # flist=filterEx(probes_expr_notf,0.5,phenoDat$`Characteristics[organism]`)
  # mirExpr_GSE32688=genes_expr_GSE32688_mir[flist,]
  
}  
  
  
# colnames(probes_expr)=substr(colnames(probes_expr),1,10)
#
#read_xls处理预后信息
phenoDat_outcome=read_xlsx("./E-MTAB-6134/E-MTAB-6134.sdrf.xlsx",sheet = 1)
phenoDat_outcome=as.data.frame(phenoDat_outcome)
#
# phenoDat=merge(phenoDat,phenoDat_outcome,by.x="patient:ch1",by.y="Tumor.ID")
phenoDat=phenoDat_outcome
rownames(phenoDat)=phenoDat$`Source Name`
## check GPL and annotate the probes to genes.
(gpl="GPL13667")
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl=gpl,type = "soft")#soft的合并性差落后，建议首选bioc
head(probe2gene)
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
# 11419 bioc
# 11571 soft 选这个
head(genes_expr)
# 修剪+median !!!!!!!!!!!!!!!!!!!!!!!!
dim(phenoDat)
phenoDat$OS.time=as.numeric(phenoDat$`Characteristics[os.delay]`)
phenoDat$OS=as.numeric(phenoDat$`Characteristics[os.event]`)
phenoDat$DFS.time=as.numeric(phenoDat$`Characteristics[dfs.delay]`)
phenoDat$DFS=as.numeric(phenoDat$`Characteristics[dfs.event]`)

# phenoDat=phenoDat[phenoDat$time>=1,]
dim(phenoDat)
# genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$`Source Name`)]
# phenoDat=phenoDat[which(phenoDat$`Source Name` %in% colnames(probes_expr)),]
genes_expr_with_clinical=genes_expr[,order(match(colnames(genes_expr),phenoDat$`Source Name`),na.last=NA)]
phenoDat=phenoDat[order(match(phenoDat$`Source Name`,colnames(genes_expr_with_clinical)),na.last=NA),]
dim(genes_expr_with_clinical)
#
library(DescTools)
genes_expr_mad_EMTAB6134=t(RobScale(t(genes_expr_with_clinical)))
genes_expr_mean_EMTAB6134=t(scale(t(genes_expr_with_clinical)))
genes_expr_EMTAB6134=genes_expr_with_clinical
probes_expr_EMTAB6134=probes_expr
phenoDat_EMTAB6134=phenoDat
save(genes_expr_EMTAB6134,genes_expr_mad_EMTAB6134,genes_expr_mean_EMTAB6134,probes_expr_EMTAB6134,phenoDat_EMTAB6134,file="EMTAB6134_after_bioc.Rdata")
write.csv(phenoDat_EMTAB6134,file="phenoDat_EMTAB6134.csv",quote = T)
pdf("E-MTAB-6134_mad.pdf",width = 100)
par(mfrow = c(2,1));
boxplot(t(genes_expr_mean_EMTAB6134))
genes_expr_mean=t(scale(t(genes_expr)))
boxplot(genes_expr_mean)
dev.off()
