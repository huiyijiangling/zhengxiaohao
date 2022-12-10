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
 
# kkk=exprs(affyRaw)
# kkk=normalizeQuantiles(kkk)
# kkk=log2(kkk+1)
# boxplot(kkk,las=2)
#gset=AnnoProbe::geoChina('GSE27342')
load("./GSE27342_eSet.Rdata")#log不平，但像是弄过了
# check the ExpressionSet
eSet=gset[[1]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)#原矩阵有错
head(probes_expr[,1:4])
probes_expr=log2(probes_expr+1)
# probes_expr=log2(probes_expr+1)
# boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])


if(eSet@annotation=="GPL570"){
system("tar -xvf ./GSE27342/GSE27342_RAW.tar -C ./GSE27342/")
celFiles <- list.celfiles('./GSE27342/',full.name=TRUE,listGzipped = T)
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
GSE27342_outcome=read.csv("./GSE27342/GSE27342_outcome.csv",header = T)
GSE27342_outcome=as.data.frame(GSE27342_outcome)
#
pdf("GSE27342_boxplot.pdf",width=100)
par(mfrow = c(2,1));
boxplot(probes_expr,las=2)
boxplot(probes_expr_rma,las=2)
dev.off()
#选择 注意每次要改 开始修剪了!!!!!!!!!
probes_expr=probes_expr#_rma
# colnames(probes_expr)=substr(colnames(probes_expr),1,10)
# phenoDat=merge(phenoDat,GSE27342_outcome,by.x="patient:ch1",by.y="Tumor.ID")
# rownames(phenoDat)=phenoDat$geo_accession
## check GPL and annotate the probes to genes.
}

#for depc
#individuals=factor(unlist(lapply(phenoDat$title,function(x) strsplit(as.character(x)," TY")[[1]][2])))
library(GDCRNATools)
# phenoDat$patientidtn=unlist(lapply(phenoDat$title,function(x) strsplit(as.character(x)," TY")[[1]][2]))
# phenoDat$patientid=substr(phenoDat$patientidtn,1,5)
phenoDat$tn=substr(phenoDat$`tissue:ch1`,1,1)

source("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/filter_f1000.R")
flist=filterEx(probes_expr,0.25,phenoDat$source_name_ch1)#??????
probes_expr=probes_expr[flist,]
library(limma)
phenoDat$patientid=unlist(lapply(phenoDat$title, function(x) rev(strsplit(x," ")[[1]])[1]))
# phenoDat$tn=ifelse(phenoDat$`tissue:ch1`=='Gastric tumor',"T","N")
probes_expr=probes_expr[,match(colnames(probes_expr),phenoDat$geo_accession)]
probes_expr=as.data.frame(probes_expr)
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
DEGAll_GSE27342 <- Limma_microarray_paired(eset     = probes_expr, 
                                           group      = phenoDat$tn, 
                                           paired     =phenoDat$patientid,
                                           comparison = 'g-n')

(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)

GPL=getGEO(gpl,AnnotGPL = F,GSEMatrix = F, getGPL =T,parseCharacteristics =T)#注意是不是下载完整了
GPL=GPL@dataTable@table

# #删繁就简 这里检验是不是注释了以前没有的新序列
# GPL[GPL==""] <- NA
# table(GPL[!is.na(GPL$miRNA_ID_LIST) & is.na(GPL$Accession),])#要有一起有
# miRNANames=GPL[is.na(GPL$Accession),]$SPOT_ID
# taa = miRNAVersionConvert(miRNANames,targetVersion = "v21",exact = TRUE)
# table(taa$TargetName)
# # 分列名称和accession号
# # 注意一定不要把ID拿去裂项，就算错了一起错还是能够配上
# GPL=separate_rows(GPL,Accession,sep = "[^-*[:alnum:].]+")
# GPL=separate_rows(GPL,miRNA_ID_LIST,sep = "[^-*[:alnum:].]+")
# GPL=unique(GPL)
# miRNANames=GPL[!is.na(GPL$miRNA_ID_LIST),]$miRNA_ID_LIST
# version=checkMiRNAVersion(miRNANames, verbose = FALSE)
# version#12
# tbb= miRNA_NameToAccession(miRNANames,version = version)#只有一个了
# tbb=unique(tbb)
# tbb2=miRNA_AccessionToName(tbb$Accession,targetVersion = "v21")
# tbb3=merge(tbb,tbb2,by="Accession",all=T)


probe2gene=GPL
probe2gene=probe2gene[!is.na(probe2gene$gene_assignment),]
probe2gene$gene_assignment=ifelse(probe2gene$gene_assignment=="---",NA,probe2gene$gene_assignment)
probe2gene=probe2gene[!is.na(probe2gene$gene_assignment),]
probe2gene$TargetName=stringr::str_split(probe2gene$gene_assignment,' // ',simplify = T)[,2]
probe2gene=unique(probe2gene[,c("ID","TargetName")])
source("./updateName.R")# ncbi与biomart相比多map了200，等于是，ncbi》biomart》org.hs
forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$TargetName)
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
# probe2gene=unique(probe2gene[,c("probe_id","ENSEMBL")])#bioc
probe2gene=unique(probe2gene[,c("ID","ENSEMBL")])#soft
# probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
colnames(probe2gene)=c("probe_id","symbol")
genes_expr <- filterEM(probes_expr,probe2gene)

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
  GSE27342_5050=TumorOnlysplitinto5050(probes_expr,"2845699",phenoDat,"geo_accession","tn","g") #ENSG00000113504
  DEGAll_GSE27342_5050 <- Limma_microarray_rma(eset     = GSE27342_5050[[1]][[2]], 
                                               group      = GSE27342_5050[[1]][[1]]$grouphl, 
                                               comparison = 'H-L')
}
DEGAll_GSE27342 <- filterEM(DEGAll_GSE27342,probe2gene)
DEGAll_GSE27342_5050 <- filterEM(DEGAll_GSE27342_5050,probe2gene)

save(DEGAll_GSE27342,DEGAll_GSE27342_5050,file ="DEGAll_GSE27342_filter.Rdata")
write.csv(DEGAll_GSE27342,quote = T,file = "DEGAll_GSE27342.csv")
write.csv(DEGAll_GSE27342_5050,quote = T,file = "DEGAll_GSE27342_5050.csv")


# 修剪+median !!!!!!!!!!!!!!!!!!!!!!!!
dim(phenoDat)
phenoDat=phenoDat[which(phenoDat$tn=="T"),]
phenoDat$name=substr(unlist(lapply(phenoDat$title,function(x) strsplit(as.character(x)," ")[[1]][6])),1,7)
phenoDat=merge(phenoDat,GSE27342_outcome,by.x="name",by.y="ID")
rownames(phenoDat)=phenoDat$geo_accession
phenoDat$time=as.numeric(phenoDat$time)
phenoDat=phenoDat[phenoDat$time>=1,]
dim(phenoDat)
genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
dim(genes_expr_with_clinical)
#
library(DescTools)
genes_expr_mad_GSE27342=t(RobScale(t(genes_expr_with_clinical)))
genes_expr_mean_GSE27342=t(scale(t(genes_expr_with_clinical)))
probes_expr_GSE27342=probes_expr
phenoDat_GSE27342=phenoDat
save(genes_expr_mad_GSE27342,genes_expr_mean_GSE27342,probes_expr_GSE27342,phenoDat_GSE27342,file="GSE27342_after_bioc.Rdata")

pdf("GSE27342_mad.pdf",width = 100)
par(mfrow = c(2,1));
boxplot(genes_expr_mad)
genes_expr_mean=t(scale(t(genes_expr)))
boxplot(genes_expr_mean)
dev.off()

DEGAll_GSE27342["ENSG00000113504",]
