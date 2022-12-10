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
# gset=AnnoProbe::geoChina('GSE26899')
load("./GSE26899_eSet.Rdata")#log不平，但像是弄过了
# check the ExpressionSet
eSet=gset[[1]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
# probes_expr=log2(probes_expr+1)

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
# phenoDat$event=substr(phenoDat$characteristics_ch1.5,7,nchar(phenoDat$characteristics_ch1.5))
# phenoDat$time=substr(phenoDat$characteristics_ch1.6,27,nchar(phenoDat$characteristics_ch1.6))
#depc
phenoDat=subset(phenoDat,phenoDat$characteristics_ch1.2!="tumor type: GIST")
Limma_microarray_rma <- function(eset, group, comparison, method='limma') {
  library(limma)
  group <- factor(group)
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  contrast.matrix <- makeContrasts(contrasts=comparison, 
                                   levels=design)
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
library(GDCRNATools)
phenoDat$tn=ifelse(is.na(phenoDat$`tumor type:ch1`),"N","T")
genes_expr=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
DEGAll_GSE26899 <- Limma_microarray_rma(eset     = genes_expr, 
                                        group      = phenoDat$tn, 
                                        comparison = 'T-N')
save(DEGAll_GSE26899,file ="DEGAll_GSE26899.Rdata" )
#############################################这里没做了，做了一半
if(eSet@annotation=="GPL6947"){

  eset <- rma(affyRaw,normalize=F)#不存在直接可以复现的quantile，因为这是用C语言写的
  # 检查数据
  probes_expr_rma=exprs(eset)
  dim(probes_expr_rma)
  #
  phenoDat_rma <- pData(eset)
  head(phenoDat_rma[,1:4])
  #
  GSE26899_outcome=read_xls("./GSE26899/GSE26899_outcome.xls",col_names = TRUE)
  GSE26899_outcome=as.data.frame(GSE26899_outcome)
  GSE26899_outcome=GSE26899_outcome[-which(substr(GSE26899_outcome$`GSM ID`,1,1)=="*"),]
  
  save(file="GSE26899_after.Rdata",probes_expr,phenoDat,probes_expr_rma,phenoDat_rma,GSE26899_outcome)
}
# 修剪+median !!!!!!!!!!!!!!!!!!!!!!!!
dim(phenoDat)
phenoDat$time=as.numeric(phenoDat$time)
phenoDat=phenoDat[phenoDat$time>=1,]
dim(phenoDat)
genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
dim(genes_expr_with_clinical)
#
library(DescTools)
genes_expr_mad_GSE26899=t(RobScale(t(genes_expr_with_clinical)))
genes_expr_mean_GSE26899=t(scale(t(genes_expr_with_clinical)))
probes_expr_GSE26899=probes_expr
phenoDat_GSE26899=phenoDat

save(genes_expr_mad_GSE26899,genes_expr_mean_GSE26899,probes_expr_GSE26899,phenoDat_GSE26899,file="GSE26899_after_bioc.Rdata")
#
if(eSet@annotation=="GPL570"){
  system("tar -xvf ./GSE26899/GSE26899_RAW.tar -C ./GSE26899/")
  celFiles <- list.celfiles('./GSE26899/',full.name=TRUE,listGzipped = T)
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
  GSE26899_outcome=read_xls("./GSE26899/GSE26899_outcome.xls",col_names = TRUE)
  GSE26899_outcome=as.data.frame(GSE26899_outcome)
  GSE26899_outcome=GSE26899_outcome[-which(substr(GSE26899_outcome$`GSM ID`,1,1)=="*"),]
  
  save(file="GSE26899_after.Rdata",probes_expr,phenoDat,probes_expr_rma,phenoDat_rma,GSE26899_outcome)
}
#
pdf("GSE26899_boxplot.pdf",width=100)
par(mfrow = c(2,1));
boxplot(probes_expr,las=2)
boxplot(probes_expr_rma,las=2)
dev.off()
