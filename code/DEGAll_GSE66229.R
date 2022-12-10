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
#gset=AnnoProbe::geoChina('GSE66229')
load("./GSE66229_eSet.Rdata")#log不平，但像是弄过了
# check the ExpressionSet
eSet=gset[[1]]#eSet=gset[[2]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)#原矩阵有错
head(probes_expr[,1:4])
# probes_expr=log2(probes_expr+1)
# boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
if(eSet@annotation=="GPL570"){
  system("tar -xvf ./GSE66229/GSE66229_RAW.tar -C ./GSE66229/")
  celFiles <- list.celfiles('./GSE66229/',full.name=TRUE,listGzipped = T)
  affyRaw <- read.celfiles(celFiles)
  # 提取矩阵并做normalization 
  #rma 自带 背景校正方法：rma，标准化方法：quantile，汇总方法：medianpolish
  #不需要重复quantile，当然limma要干嘛随便
  eset <- rma(affyRaw)#不存在直接可以复现的quantile，因为这是用C语言写的
  # 检查数据
  probes_expr_rma=exprs(eset)
  dim(probes_expr_rma)
  #
  pdf("GSE66229_boxplot.pdf",width=100)
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
# phenoDat_outcome=read.csv("./GSE66229/GSE66229_outcome.csv",header = T)
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
# 发现有可能是一对多
probe2gene=separate_rows(probe2gene,symbol,sep = " /// ")
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

Limma_microarray_paired <- function(eset, group, comparison,paired,method='limma') {
  library(limma)#一应注意limma 的makecontrast和level选一个，否则顺序按字母,分组是必须是对照组在前，实验组在后
  group <- factor(group,levels = rev(strsplit(comparison,"-")[[1]]),ordered = F)
  paired <- factor(paired) 
  design <- model.matrix(~paired+group)
  fit <- lmFit(eset, design)
  fit2 <- eBayes(fit)
  DEGAll <- topTable(fit2, coef=paste0("group",strsplit(comparison,"-")[[1]][1]), n = Inf)
  colnames(DEGAll) <- c('logFC', 'AveExpr', 't', 
                        'PValue', 'FDR', 'B')
  DEGAll$FDR <- p.adjust(DEGAll$PValue, method = 'fdr')
  o <- order(DEGAll$FDR)
  DEGAll <- DEGAll[o,]
  
  if (startsWith(rownames(eset)[1], 'ENSG')) {
    degList <- biotype[match(rownames(DEGAll), biotype$ensemblID),]
    degOutput <- data.frame(symbol=degList$geneSymbol, 
                            group=degList$group, DEGAll)
    
    keep <- which(! is.na(degOutput$symbol))
    degOutput <- degOutput[keep,]
    return(degOutput)
  } else {
    return (DEGAll)
  }
}
Limma_microarray_rma <- function(eset, group, comparison, method='limma') {
  library(limma)
  group <- factor(group)
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  contrast.matrix <- makeContrasts(contrasts=comparison, 
                                   levels=design)#将实验组-对照组。coef选1，可无视字母顺序
  # contrast.matrix<-makeContrasts(paste0(unique(group_list),
  #                                       collapse = "-"),levels = design)
  fit <- lmFit(eset, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  DEGAll <- topTable(fit2, coef=1, n = Inf)
  colnames(DEGAll) <- c('logFC', 'AveExpr', 't', 
                        'PValue', 'FDR', 'B')
  
  DEGAll$FDR <- p.adjust(DEGAll$PValue, method = 'fdr')
  o <- order(DEGAll$FDR)
  DEGAll <- DEGAll[o,]
  
  # if (startsWith(rownames(eset)[1], 'ENSG')) {
  #   degList <- biotype[match(rownames(DEGAll), biotype$ensemblID),]
  #   degOutput <- data.frame(symbol=degList$geneSymbol, 
  #                           group=degList$group, DEGAll)
  #   
  #   keep <- which(! is.na(degOutput$symbol))
  #   degOutput <- degOutput[keep,]
  #   return(degOutput)
  # } else {
  return (DEGAll)
  # }
}
library(GDCRNATools)
colnames(genes_expr)=stringr::str_split(colnames(genes_expr),"_",simplify = T)[,1]
phenoDat=phenoDat[match(colnames(genes_expr),phenoDat$geo_accession),]
phenoDat$tn=ifelse(phenoDat$`tissue:ch1`=='Gastric tumor',"T","N")
#b[match(A$id,B$id),]#从B(10)里面提取A(2)，按照A的顺序,想用match排序的化注意
DEGAll_GSE66229 <- Limma_microarray_rma(eset     = genes_expr, 
                                        group      = phenoDat$tn, 
                                        comparison = 'T-N')
# wyz n01
load("phenoDat_GSE62254_wyzpn01.Rdata")
phenoDat_GSE62254$pn01=ifelse(phenoDat_GSE62254$pn01==1,"P","N")
genes_expr_T=subset(genes_expr,select=phenoDat_GSE62254$geo_accession)
DEN01_GSE66229 <- Limma_microarray_rma(eset     = genes_expr_T, 
                                        group      = phenoDat_GSE62254$pn01, 
                                        comparison = 'P-N')
genes_expr_T_GSE62254=genes_expr_T
save(DEN01_GSE66229,genes_expr_T_GSE62254,phenoDat_GSE62254,file = "DEN01_GSE66229.Rdata")
if(T){
  TumorOnlysplitinto5050 <- function(datasets,selectgroup,phenoDat,coln,traitV,trait){
    new_phenoDat=subset(phenoDat,phenoDat[[traitV]]==trait)
    new_datasets=as.data.frame(datasets)
    new_datasets=subset(new_datasets,select=colnames(new_datasets) %in% new_phenoDat[[coln]])
    try(if(anyNA(new_datasets[selectgroup,])) stop("Please impute NA value first!"))
    # new_datasets$grouphl=ifelse(new_datasets[[selectgroup]]>quantile(new_datasets[[selectgroup]])[3],1,
    #                             ifelse(new_datasets[[selectgroup]]<=quantile(new_datasets[[selectgroup]])[3],0,NA))
    # df[,-which(names(df)%in%c("a","b")]
    new_datasetsA=subset(new_datasets,select=new_datasets[selectgroup,]>(quantile(new_datasets[selectgroup,])[[3]]))
    new_phenoDatA=new_phenoDat[match(colnames(new_datasetsA),new_phenoDat[[coln]]),]
    new_datasetsB=subset(new_datasets,select=new_datasets[selectgroup,]<=(quantile(new_datasets[selectgroup,])[[3]]))
    new_phenoDatB=new_phenoDat[match(colnames(new_datasetsB),new_phenoDat[[coln]]),]
    new_datasets=cbind(new_datasetsA,new_datasetsB)
    new_phenoDat=rbind(new_phenoDatA,new_phenoDatB)
    new_phenoDat$grouphl=c(rep("H",nrow(new_phenoDatA)),rep("L",nrow(new_phenoDatB)))
    return(list(list(new_phenoDat,new_datasets),list(new_phenoDatA,new_datasetsA),list(new_phenoDatB,new_datasetsB)))
  }
  GSE66229_5050=TumorOnlysplitinto5050(genes_expr,"ENSG00000113504",phenoDat,"geo_accession","tn","T")
  DEGAll_GSE66229_5050 <- Limma_microarray_rma(eset     = GSE66229_5050[[1]][[2]], 
                                          group      = GSE66229_5050[[1]][[1]]$grouphl, 
                                          comparison = 'H-L')
  
  }


save(DEGAll_GSE66229,DEGAll_GSE66229_5050,file ="DEGAll_GSE66229_filter.Rdata")
write.csv(DEGAll_GSE66229,quote = T,file = "DEGAll_GSE66229.csv")
write.csv(DEGAll_GSE66229_5050,quote = T,file = "DEGAll_GSE66229_5050.csv")
# save(DEGAll_GSE66229,file ="DEGAll_GSE66229.Rdata" )

# load("DEGAll_GSE66229.Rdata")
