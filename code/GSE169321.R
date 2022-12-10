#阴性结果
#默认的panel 上传的矩阵，已经z-score了,实在不行直接scale后算了
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
# gset=AnnoProbe::geoChina('GSE169321')
# gset=getGEO('GSE169321',getGPL = F)
# save(gset,file="GSE169321_eSet.Rdata")
load("./GSE169321_eSet.Rdata")#log不平，但像是弄过了
# check the ExpressionSet
eSet=gset[[1]]#eSet=gset[[2]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)#原矩阵有错
head(probes_expr[,1:4])
probes_expr=read.csv("C:\\Users\\zxh\\Desktop\\R\\20201025 pancreatic cancer validation\\GSE169321\\GSE169321_raw_count_all.csv.gz")
rownames(probes_expr)=probes_expr$id
probes_expr$id=NULL
# probes_expr=as.data.frame(t(probes_expr))
# probes_expr=log2(probes_expr+1)
# boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
write.csv(phenoDat,file="phenoDat_GSE169321.csv",quote = T)
phenoDat$sample_type=unlist(lapply(phenoDat$title,function(x) strsplit(as.character(x)," ")[[1]][[1]]))
phenoDat=subset(phenoDat,sample_type %in% "Primary")
phenoDat$patient=paste0("T",unlist(lapply(phenoDat$title,function(x) rev(strsplit(as.character(x)," ")[[1]])[[1]])))
phenoDat[6,"patient"] <- "Primary_44T9"
probes_expr=subset(probes_expr,select=phenoDat$patient)
phenoDat$group=ifelse(phenoDat$`treatment:ch1`=="FOLFIRINOX","neo","none")
library(GDCRNATools)


probe2gene=data.frame(probe_id=rownames(probes_expr),symbol=rownames(probes_expr))
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
probes_expr=probes_expr
genes_expr <- filterEM(probes_expr,probe2gene)
genes_expr_GSE169321_rna=genes_expr





DEGAll <- gdcDEAnalysis(counts     = genes_expr_GSE169321_rna,
                        group      = phenoDat$group,
                        comparison = 'neo-none',
                        method     = 'DESeq2',
                        filter =T)#默认是true，不筛选则选择false
DEGAll=subset(DEGAll,DEGAll$group=="protein_coding")
DEGAll$FDR=as.numeric(DEGAll$FDR)
DEGAll=subset(DEGAll,DEGAll$FDR<0.05)
dakala=DEGAll$symbol
save(dakala,"dakala.Rdata")
# boxplot(log2(genes_expr_GSE169321_rna+1))
#非常好顺序不重要
#DEGAll["FCGBP",]
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all', fc = 1, pval = 0.05)#
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 1.000, pval = 0.05)#
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 1.000, pval = 0.05)#其实放到0.05也可以
LNC_0_F <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 100)#
PC_0_F <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 100)
