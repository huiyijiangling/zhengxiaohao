rm(list=ls())
options(stringsAsFactors = F)
gc()
#cols <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(2)
#注意unlist mirlist 因为里面放了/!!!!!!!!!!!!
#我们可以过滤但是不要超过50%
source("./scripts/filter_f1000.R")
library(stringr)
# load("starbase_homo_sig_cerna_all.Rdata")
library(GDCRNATools)
#use_same_name
load("PAAD_GDCRNATOOLS.Rdata")
#TCGA-HZ-A9TJ-06A	是一个meta标本已经删除

#TCGA-YH-A8SY-01 是NA

metaMatrix.RNA=subset(metaMatrix.RNA,!sample %in%c("TCGA-3A-A9IV-01","TCGA-3A-A9IN-01","TCGA-3A-A9IL-01","TCGA-3A-A9IL-01","TCGA-2L-AAQM-01","TCGA-3A-A9IO-01","TCGA-3A-A9IR-01","TCGA-3A-A9IJ-01","TCGA-3A-A9IS-01","TCGA-YH-A8SY-01","TCGA-H8-A6C1-01"))
metaMatrix.MIR=subset(metaMatrix.MIR,!sample %in%c("TCGA-3A-A9IV-01","TCGA-3A-A9IN-01","TCGA-3A-A9IL-01","TCGA-3A-A9IL-01","TCGA-2L-AAQM-01","TCGA-3A-A9IO-01","TCGA-3A-A9IR-01","TCGA-3A-A9IJ-01","TCGA-3A-A9IS-01","TCGA-YH-A8SY-01","TCGA-H8-A6C1-01"))
mirCounts=subset(mirCounts,select=metaMatrix.MIR$sample)#172=168+4
table(substr(colnames(mirCounts),14,15))
rnaCounts=subset(rnaCounts,select=metaMatrix.RNA$sample)#171=167+4
table(substr(colnames(rnaCounts),14,15))
#"TCGA-F2-6880-01","TCGA-H6-A45N-11","TCGA-2J-AABP-01"
mirCounts[,"TCGA-IB-7647-01"]
rnaCounts[,"TCGA-IB-7647-01"]
metaMatrix.MIR=metaMatrix.MIR[metaMatrix.MIR$sample %in% metaMatrix.RNA$sample,]#TCGA-IB-7647-01	是一个没有标本
metaMatrix.RNA=metaMatrix.RNA[metaMatrix.RNA$sample %in% metaMatrix.MIR$sample,]
mirCounts=subset(mirCounts,select=metaMatrix.MIR$sample)#167+4
table(substr(colnames(mirCounts),14,15))
rnaCounts=subset(rnaCounts,select=metaMatrix.RNA$sample)#167+4
table(substr(colnames(rnaCounts),14,15))

mirExpr_tcga <- gdcVoomNormalization_modify(counts = mirCounts)#limma filter
rnaExpr_tcga <- gdcVoomNormalization_modify(counts = rnaCounts)#limma filter


if(F){
  rnaCounts_rownames <- rownames(rnaCounts)
  rnaCounts_colnames <- colnames(rnaCounts)
  rnaCounts_quant <- preprocessCore::normalize.quantiles(
    as.matrix(rnaCounts))
  rownames(rnaCounts_quant) <- rnaCounts_rownames
  colnames(rnaCounts_quant) <- rnaCounts_colnames
  rnaCounts_quant=as.data.frame(rnaCounts_quant)
  rnaExpr_quant=log2(rnaCounts_quant+1)
  
  
  # rnaCounts_quant=limma::normalizeBetweenArrays(log2(rnaCounts+1))
  mirCounts_rownames <- rownames(mirCounts)
  mirCounts_colnames <- colnames(mirCounts)
  mirCounts_quant <- preprocessCore::normalize.quantiles(
    as.matrix(mirCounts))
  rownames(mirCounts_quant) <- mirCounts_rownames
  colnames(mirCounts_quant) <- mirCounts_colnames
  mirCounts_quant=as.data.frame(mirCounts_quant)
  mirExpr_quant=log2(mirCounts_quant+1) 
  
  
  mirExpr_tcga=mirExpr_quant[rownames(mirExpr_tcga),]
  rnaExpr_tcga=rnaExpr_quant[rownames(rnaExpr_tcga),]
  
}


#
load("GSE119794_before.Rdata")
mirExpr_GSE119794 <- gdcVoomNormalization_modify(counts = genes_expr_GSE119794_mir)#limma filter
rnaExpr_GSE119794<- gdcVoomNormalization_modify(counts = genes_expr_GSE119794_rna)#limma filter
if(F){
  genes_expr_GSE119794_rna_rownames <- rownames(genes_expr_GSE119794_rna)
  genes_expr_GSE119794_rna_colnames <- colnames(genes_expr_GSE119794_rna)
  genes_expr_GSE119794_rna_quant <- preprocessCore::normalize.quantiles(
    as.matrix(genes_expr_GSE119794_rna))
  rownames(genes_expr_GSE119794_rna_quant) <- genes_expr_GSE119794_rna_rownames
  colnames(genes_expr_GSE119794_rna_quant) <- genes_expr_GSE119794_rna_colnames
  genes_expr_GSE119794_rna_quant=as.data.frame(genes_expr_GSE119794_rna_quant)
  rnaExpr_quant=log2(genes_expr_GSE119794_rna_quant+1)
  
  genes_expr_GSE119794_mir_rownames <- rownames(genes_expr_GSE119794_mir)
  genes_expr_GSE119794_mir_colnames <- colnames(genes_expr_GSE119794_mir)
  genes_expr_GSE119794_mir_quant <- preprocessCore::normalize.quantiles(
    as.matrix(genes_expr_GSE119794_mir))
  rownames(genes_expr_GSE119794_mir_quant) <- genes_expr_GSE119794_mir_rownames
  colnames(genes_expr_GSE119794_mir_quant) <- genes_expr_GSE119794_mir_colnames
  genes_expr_GSE119794_mir_quant=as.data.frame(genes_expr_GSE119794_mir_quant)
  mirExpr_quant=log2(genes_expr_GSE119794_mir_quant+1) 
  
  
  mirExpr_GSE119794=mirExpr_quant[rownames(mirExpr_GSE119794),]
  rnaExpr_GSE119794=rnaExpr_quant[rownames(rnaExpr_GSE119794),]
  
}
phenoDat_GSE119794$newname=paste(phenoDat_GSE119794$`subject id:ch1`,phenoDat_GSE119794$source_name_ch1,sep = "_")
colnames(rnaExpr_GSE119794)=phenoDat_GSE119794[match(colnames(rnaExpr_GSE119794), phenoDat_GSE119794$geo_accession),c("newname")]
colnames(mirExpr_GSE119794)=phenoDat_GSE119794[match(colnames(mirExpr_GSE119794), phenoDat_GSE119794$geo_accession),c("newname")]
phenoDat_GSE119794_rna=subset(phenoDat_GSE119794,str_split(phenoDat_GSE119794$title,"_",simplify = T)[,2]=="mRNA")
phenoDat_GSE119794_mir=subset(phenoDat_GSE119794,str_split(phenoDat_GSE119794$title,"_",simplify = T)[,2]=="microRNA")
load("GSE32688_after_bioc.Rdata")#
flist=filterEx(genes_expr_GSE32688_mir,0.35,phenoDat_GSE32688_mir$source_name_ch1)
mirExpr_GSE32688=genes_expr_GSE32688_mir[flist,]
flist=filterEx(genes_expr_GSE32688_rna,0.35,phenoDat_GSE32688_rna$source_name_ch1)
rnaExpr_GSE32688=genes_expr_GSE32688_rna[flist,]
phenoDat_GSE32688_rna$newname=gsub("[^[:alnum:]]","",phenoDat_GSE32688_rna$title)
phenoDat_GSE32688_mir$newname=gsub("[^[:alnum:]]","",phenoDat_GSE32688_mir$title)
colnames(rnaExpr_GSE32688)=phenoDat_GSE32688_rna[match(colnames(rnaExpr_GSE32688), phenoDat_GSE32688_rna$geo_accession),c("newname")]
colnames(mirExpr_GSE32688)=phenoDat_GSE32688_mir[match(colnames(mirExpr_GSE32688), phenoDat_GSE32688_mir$geo_accession),c("newname")]
load("GSE43797_after_bioc.Rdata")
flist=filterEx(genes_expr_GSE43797_mir,0.5,phenoDat_GSE43797_mir$source_name_ch1)
mirExpr_GSE43797=genes_expr_GSE43797_mir[flist,]
flist=filterEx(genes_expr_GSE43797_rna,0.5,phenoDat_GSE43797_rna$source_name_ch1)
rnaExpr_GSE43797=genes_expr_GSE43797_rna[flist,]
phenoDat_GSE43797_rna$newname=gsub("[^[:alnum:]]","",phenoDat_GSE43797_rna$title)
phenoDat_GSE43797_mir$newname=str_split(phenoDat_GSE43797_mir$title," ",simplify = T)[,1]
phenoDat_GSE43797_mir$newname=gsub("[^[:alnum:].]+","",phenoDat_GSE43797_mir$newname)
colnames(rnaExpr_GSE43797)=phenoDat_GSE43797_rna[match(colnames(rnaExpr_GSE43797), phenoDat_GSE43797_rna$geo_accession),c("newname")]
colnames(mirExpr_GSE43797)=phenoDat_GSE43797_mir[match(colnames(mirExpr_GSE43797), phenoDat_GSE43797_mir$geo_accession),c("newname")]
# load("GSE41372_after_bioc.Rdata")#唯独这个集合用bioc 更少了
load("GSE41372_after_soft.Rdata")
# boxplot(log2(probes_expr_GSE41372_rna+1))#62-64 log2后6
# boxplot(log2(probes_expr_GSE41372_mir+1))#19 log后4
genes_expr_GSE41372_mir=log2(genes_expr_GSE41372_mir+1)
genes_expr_GSE41372_rna=log2(genes_expr_GSE41372_rna+1)
flist=filterEx(genes_expr_GSE41372_mir,0.25,phenoDat_GSE41372_mir$`tissue:ch1`)
mirExpr_GSE41372=genes_expr_GSE41372_mir[flist,]
flist=filterEx(genes_expr_GSE41372_rna,0.25,phenoDat_GSE41372_rna$`tissue:ch1`)
rnaExpr_GSE41372=genes_expr_GSE41372_rna[flist,]
phenoDat_GSE41372_rna$newname=gsub("[^[:alnum:]]","",phenoDat_GSE41372_rna$title)
phenoDat_GSE41372_mir$newname=gsub("[^[:alnum:]]","",phenoDat_GSE41372_mir$title)
colnames(rnaExpr_GSE41372)=phenoDat_GSE41372_rna[match(colnames(rnaExpr_GSE41372), phenoDat_GSE41372_rna$geo_accession),c("newname")]
colnames(mirExpr_GSE41372)=phenoDat_GSE41372_mir[match(colnames(mirExpr_GSE41372), phenoDat_GSE41372_mir$geo_accession),c("newname")]
#其实需要看一看是不是已经把低表达单位filter掉了
if(F){
  library(miRBaseConverter)
  miRNANames=rownames(as.data.frame(mirCounts))
  version=checkMiRNAVersion(miRNANames, verbose = FALSE)
  version#21
}
#+scale_y_continuous(breaks=seq(from=-4,to=16)) #seq(-2,0,2,4,6,8,10,12,14,16) ,breaks=seq(from=-4,to=16,by=1)
if(F){
  pdf("GSE32688_mad.pdf",width = 100)
  par(mfrow = c(4,1));
  boxplot(probes_expr_GSE32688_rna)#3
  boxplot(probes_expr_GSE32688_mir)#5
  # genes_expr_mean=t(scale(t(genes_expr)))
  # boxplot(genes_expr_mean)
  boxplot(genes_expr_mean_GSE32688_rna)
  boxplot(genes_expr_mean_GSE32688_mir)
  dev.off()
  
  pdf("GSE43797_mad.pdf",width = 100)
  par(mfrow = c(4,1));
  boxplot(probes_expr_GSE43797_rna)#6.8 gaole
  boxplot(probes_expr_GSE43797_mir)#5.6
  # genes_expr_mean=t(scale(t(genes_expr)))
  # boxplot(genes_expr_mean)
  boxplot(genes_expr_mean_GSE43797_rna)
  boxplot(genes_expr_mean_GSE43797_mir)
  dev.off()
  
  pdf("GSE41372_mad.pdf",width = 100)
  par(mfrow = c(4,1));
  boxplot(probes_expr_GSE41372_rna)#62-64 log2后6
  boxplot(probes_expr_GSE41372_mir)#19 log后4
  # genes_expr_mean=t(scale(t(genes_expr)))
  # boxplot(genes_expr_mean)
  boxplot(genes_expr_mean_GSE41372_rna)
  boxplot(genes_expr_mean_GSE41372_mir)
  dev.off()
  
  pdf("GSE119794_mad.pdf",width = 100)
  par(mfrow = c(4,1));
  boxplot(probes_expr_GSE119794_rna)# 80 log 6  1/2/4表达很低
  boxplot(probes_expr_GSE119794_mir)# 01??? 0001
  # genes_expr_mean=t(scale(t(genes_expr)))
  # boxplot(genes_expr_mean)
  boxplot(genes_expr_mean_GSE119794_rna)
  boxplot(genes_expr_mean_GSE119794_mir)
  dev.off()
}

common_ensg=list(TCGA=rownames(as.data.frame(rnaExpr_tcga)),
                 GSE32688=rownames(rnaExpr_GSE32688),
                 GSE119794=rownames(rnaExpr_GSE119794),
                 GSE41372=rownames(rnaExpr_GSE41372),
                 GSE43797=rownames(rnaExpr_GSE43797)
)
common_ensg=Reduce(intersect,common_ensg)

if(T){
  library(tidyr)
  seprowname<- function(datasets){
    datasets=as.data.frame(datasets)
    datasets$rn=rownames(as.data.frame(datasets))
    datasets=separate_rows(datasets,rn,sep = "/")
    datasets=as.data.frame(datasets)
    rownames(datasets)=datasets$rn
    datasets=datasets[,-which(colnames(datasets) %in% c("rn"))]
    return(datasets)
  }
  mirExpr_tcga=seprowname(mirExpr_tcga)
  mirExpr_GSE32688=seprowname(mirExpr_GSE32688)
  mirExpr_GSE119794=seprowname(mirExpr_GSE119794)
  mirExpr_GSE41372=seprowname(mirExpr_GSE41372)
  mirExpr_GSE43797=seprowname(mirExpr_GSE43797)
}
common_mir=list(TCGA=rownames(as.data.frame(mirExpr_tcga)),
                GSE32688=rownames(mirExpr_GSE32688),
                GSE119794=rownames(mirExpr_GSE119794),
                GSE41372=rownames(mirExpr_GSE41372),
                GSE43797=rownames(mirExpr_GSE43797)
)
common_mir=Reduce(intersect,common_mir) 
#
load("genecodev36.Rdata")
common_ensg_v36=genecodev36[genecodev36$gene_id %in% common_ensg,]
table(common_ensg_v36$gene_type)
common_ensg_v36=common_ensg_v36[common_ensg_v36$gene_type %in% c('lncRNA','protein_coding'),]
mrna=common_ensg_v36[common_ensg_v36$gene_type %in% c("protein_coding"),]
lncrna=common_ensg_v36[common_ensg_v36$gene_type %in% c("lncRNA"),]
load(file = "starbase_homo_sig_cerna_all.Rdata")
rnainter_homo_sig_cernalnc=rnainter_homo_sig_cernalnc[rnainter_homo_sig_cernalnc$miRNAname %in% common_mir,]
rnainter_homo_sig_cernamrna=rnainter_homo_sig_cernamrna[rnainter_homo_sig_cernamrna$miRNAname %in% common_mir,]
rnainter_homo_sig_cernalnc=rnainter_homo_sig_cernalnc[rnainter_homo_sig_cernalnc$geneID %in% common_ensg_v36$gene_id,]
rnainter_homo_sig_cernamrna=rnainter_homo_sig_cernamrna[rnainter_homo_sig_cernamrna$geneID %in% common_ensg_v36$gene_id,]
lnc_inter_target=split(rnainter_homo_sig_cernalnc[,1], as.character(rnainter_homo_sig_cernalnc$geneID))
mrna_inter_target=split(rnainter_homo_sig_cernamrna[,1], as.character(rnainter_homo_sig_cernamrna$geneID))




if(F){
  ceOutput_tcga_all <- gdcCEAnalysis(lnc         = lncrna$gene_id, 
                                     pc          = mrna$gene_id,
                                     deMIR       = common_mir,
                                     lnc.targets = lnc_inter_target,#"starBase",#lnc_inter_target
                                     pc.targets  = mrna_inter_target,#"starBase",#mrna_inter_target
                                     rna.expr    = rnaExpr_tcga, 
                                     mir.expr    = mirExpr_tcga)#超几何分布及p值
  ceOutput_GSE119794_all <- gdcCEAnalysis(lnc         = lncrna$gene_id, 
                                          pc          = mrna$gene_id,
                                          deMIR       = common_mir,
                                          lnc.targets = lnc_inter_target,#"starBase",#lnc_inter_target
                                          pc.targets  = mrna_inter_target,#"starBase",#mrna_inter_target
                                          rna.expr    = rnaExpr_GSE119794, 
                                          mir.expr    = mirExpr_GSE119794)#超几何分布及p值
  ceOutput_GSE43797_all <- gdcCEAnalysis(lnc         = lncrna$gene_id, 
                                         pc          = mrna$gene_id,
                                         deMIR       = common_mir,
                                         lnc.targets = lnc_inter_target,#"starBase",#lnc_inter_target
                                         pc.targets  = mrna_inter_target,#"starBase",#mrna_inter_target
                                         rna.expr    = rnaExpr_GSE43797, 
                                         mir.expr    = mirExpr_GSE43797)#超几何分布及p值
  ceOutput_GSE41372_all <- gdcCEAnalysis(lnc         = lncrna$gene_id, 
                                         pc          = mrna$gene_id,
                                         deMIR       = common_mir,
                                         lnc.targets = lnc_inter_target,#"starBase",#lnc_inter_target
                                         pc.targets  = mrna_inter_target,#"starBase",#mrna_inter_target
                                         rna.expr    = rnaExpr_GSE41372, 
                                         mir.expr    = mirExpr_GSE41372)#超几何分布及p值
  ceOutput_GSE32688_all <- gdcCEAnalysis(lnc         = lncrna$gene_id, 
                                         pc          = mrna$gene_id,
                                         deMIR       = common_mir,
                                         lnc.targets = lnc_inter_target,#"starBase",#lnc_inter_target
                                         pc.targets  = mrna_inter_target,#"starBase",#mrna_inter_target
                                         rna.expr    = rnaExpr_GSE32688, 
                                         mir.expr    = mirExpr_GSE32688)#超几何分布及p值
  
  ceOutput_tcga_all$preserved=paste(ceOutput_tcga_all$lncRNAs,ceOutput_tcga_all$Genes,sep="_")
  ceOutput_GSE119794_all$preserved=paste(ceOutput_GSE119794_all$lncRNAs,ceOutput_GSE119794_all$Genes,sep="_")
  ceOutput_GSE43797_all$preserved=paste(ceOutput_GSE43797_all$lncRNAs,ceOutput_GSE43797_all$Genes,sep="_")
  ceOutput_GSE41372_all$preserved=paste(ceOutput_GSE41372_all$lncRNAs,ceOutput_GSE41372_all$Genes,sep="_")
  ceOutput_GSE32688_all$preserved=paste(ceOutput_GSE32688_all$lncRNAs,ceOutput_GSE32688_all$Genes,sep="_")
  
  ceOutput_tcga_all_co=subset(ceOutput_tcga_all,hyperPValue<0.05&cor>=0.0&corPValue<0.05)
  ceOutput_GSE119794_all_co=subset(ceOutput_GSE119794_all,hyperPValue<0.05&cor>=0.0&corPValue<0.05)
  ceOutput_GSE43797_all_co=subset(ceOutput_GSE43797_all,hyperPValue<0.05&cor>=0.0&corPValue<0.05)
  ceOutput_GSE41372_all_co=subset(ceOutput_GSE41372_all,hyperPValue<0.05&cor>=0.0&corPValue<0.05)
  ceOutput_GSE32688_all_co =subset(ceOutput_GSE32688_all,hyperPValue<0.05&cor>=0.0&corPValue<0.05)
  
  ceOutput_all_co=list(TCGA=ceOutput_tcga_all_co$preserved,
                       GSE119794=ceOutput_GSE119794_all_co$preserved,
                       GSE43797=ceOutput_GSE43797_all_co$preserved,
                       GSE41372=ceOutput_GSE41372_all_co$preserved,
                       GSE32688=ceOutput_GSE32688_all_co$preserved)
  ceOutput_all_co=Reduce(intersect,ceOutput_all_co) 
  ceOutput_all_co=as.data.frame(ceOutput_all_co)
  fnmdx25=subset(ceOutput_tcga_all_co,preserved %in% ceOutput_all_co$ceOutput_all_co)
  write.csv(fnmdx25,file = "tab.csv",quote = T)
}


#
all_ensg=cbind(rnaExpr_tcga[common_ensg,],
               rnaExpr_GSE32688[common_ensg,],
               rnaExpr_GSE119794[common_ensg,],
               rnaExpr_GSE41372[common_ensg,],
               rnaExpr_GSE43797[common_ensg,])

all_ensg_phe=rbind(as.matrix(metaMatrix.RNA[,c("sample","sample_type")]),
                   as.matrix(phenoDat_GSE32688_rna[,c("newname","disease status:ch1")]),
                   as.matrix(phenoDat_GSE119794_rna[,c("newname","source_name_ch1")]),
                   as.matrix(phenoDat_GSE41372_rna[,c("newname","tissue:ch1")]),
                   as.matrix(phenoDat_GSE43797_rna[,c("newname","tissue:ch1")]))
all_ensg_phe=as.data.frame(all_ensg_phe)
all_ensg_phe=cbind(all_ensg_phe,datasets = c(rep("TCGA",nrow(metaMatrix.RNA)),
                                             rep("GSE32688",nrow(phenoDat_GSE32688_rna)),
                                             rep("GSE119794",nrow(phenoDat_GSE119794_rna)),
                                             rep("GSE41372",nrow(phenoDat_GSE41372_rna)),
                                             rep("GSE43797",nrow(phenoDat_GSE43797_rna))
))


all_ensg_phe$TN=ifelse(all_ensg_phe$sample_type %in% c("PrimaryTumor","pancreatic cancer","PC","pancreatic ductal adenocarcinoma","ductal adenocarcinoma (PCA)"),"Tumor",ifelse(all_ensg_phe$sample_type %in% c("SolidTissueNormal","non-malignant pancreas","Normal","normal pancreas","non-neoplastic pancreas"),"Normal","wrong"))
all_ensg_phe$subset=paste0(all_ensg_phe$datasets,"_",all_ensg_phe$TN)
all_ensg=all_ensg[,all_ensg_phe$sample]                      

table(all_ensg_phe$sample_type)

#as.matrix可以rbind
# 名字同原来已有的名字不相对
##示例
if(!require(factoextra))install.packages("factoextra")
if(!require(FactoMineR))install.packages("FactoMineR")
library(factoextra)
library(FactoMineR)
# "all" "none"

probesetvar = apply(all_ensg, 1, var)  #表达变化大的基因
ord = order(probesetvar, decreasing=TRUE)[1:500]   #前200个基因，或者更多
pca = prcomp(t(all_ensg[ord,]), scale=TRUE)
ss=summary(pca)
# 样本聚类

pdf(file='pca_before_TN.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_ensg_phe$TN),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1) 
  #legend.position = "top",legend.direction = "horizontal",legend.box = "vertical",
  # scale_color_manual(values = cols) +
  # scale_fill_manual(values = cols) +
dev.off()
pdf(file='pca_before_series.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_ensg_phe$subset),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='pca_before_studies.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_ensg_phe$datasets),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
# 4sigma 0.999936657516334
# 6sigma 0.999999998026825

# ,mod = NULL,这里的含义是协变量，告诉其本身有内在的分组，如果不删除tcga里癌旁正常，可以留着。 一定是先校正后删除样本 20201020
library(sva)
mod = model.matrix(~as.factor(TN), data=all_ensg_phe)
all_ensg=ComBat(dat=all_ensg,batch=all_ensg_phe$datasets,mod=mod)
# rnaExpr=as.data.frame(rnaExpr_quant) 

probesetvar = apply(all_ensg, 1, var)  #表达变化大的基因
ord = order(probesetvar, decreasing=TRUE)[1:500]   #前200个基因，或者更多
pca = prcomp(t(all_ensg[ord,]), scale=TRUE)
ss=summary(pca)
pdf(file='pca_after_TN.pdf',height = 8,width = 8)#var,quali,all,none
fviz_pca_ind(pca, label="none", habillage=as.factor(all_ensg_phe$TN),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='pca_after_series.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_ensg_phe$subset),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='pca_after_studies.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_ensg_phe$datasets),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()


#"TCGA-F2-6880-01","TCGA-H6-A45N-11","TCGA-2J-AABP-01"
metaMatrix.RNA=subset(metaMatrix.RNA,!sample %in%c("TCGA-F2-6880-01","TCGA-H6-A45N-11","TCGA-2J-AABP-01"))
rnaExpr_tcga=rnaExpr_tcga[,metaMatrix.RNA$sample]
phenoDat_GSE32688_rna=subset(phenoDat_GSE32688_rna,!newname %in%c("humanPDAC24"))
rnaExpr_GSE32688=rnaExpr_GSE32688[,phenoDat_GSE32688_rna$newname]
#
all_ensg=cbind(rnaExpr_tcga[common_ensg,],
               rnaExpr_GSE32688[common_ensg,],
               rnaExpr_GSE119794[common_ensg,],
               rnaExpr_GSE41372[common_ensg,],
               rnaExpr_GSE43797[common_ensg,])

all_ensg_phe=rbind(as.matrix(metaMatrix.RNA[,c("sample","sample_type")]),
                   as.matrix(phenoDat_GSE32688_rna[,c("newname","disease status:ch1")]),
                   as.matrix(phenoDat_GSE119794_rna[,c("newname","source_name_ch1")]),
                   as.matrix(phenoDat_GSE41372_rna[,c("newname","tissue:ch1")]),
                   as.matrix(phenoDat_GSE43797_rna[,c("newname","tissue:ch1")]))
all_ensg_phe=as.data.frame(all_ensg_phe)
all_ensg_phe=cbind(all_ensg_phe,datasets = c(rep("TCGA",nrow(metaMatrix.RNA)),
                                             rep("GSE32688",nrow(phenoDat_GSE32688_rna)),
                                             rep("GSE119794",nrow(phenoDat_GSE119794_rna)),
                                             rep("GSE41372",nrow(phenoDat_GSE41372_rna)),
                                             rep("GSE43797",nrow(phenoDat_GSE43797_rna))
))
all_ensg_phe$TN=ifelse(all_ensg_phe$sample_type %in% c("PrimaryTumor","pancreatic cancer","PC","pancreatic ductal adenocarcinoma","ductal adenocarcinoma (PCA)"),"Tumor",ifelse(all_ensg_phe$sample_type %in% c("SolidTissueNormal","non-malignant pancreas","Normal","normal pancreas","non-neoplastic pancreas"),"Normal","wrong"))
all_ensg_phe$subset=paste0(all_ensg_phe$datasets,"_",all_ensg_phe$TN)
all_ensg=all_ensg[,all_ensg_phe$sample]                      
library(sva)
mod = model.matrix(~as.factor(TN), data=all_ensg_phe)
all_ensg=ComBat(dat=all_ensg,batch=all_ensg_phe$datasets,mod=mod)
# rnaExpr=as.data.frame(rnaExpr_quant) 

probesetvar = apply(all_ensg, 1, var)  #表达变化大的基因
ord = order(probesetvar, decreasing=TRUE)[1:500]   #前200个基因，或者更多
pca = prcomp(t(all_ensg[ord,]), scale=TRUE)
ss=summary(pca)
pdf(file='final_pca_after_TN.pdf',height = 8,width = 8)#var,quali,all,none
fviz_pca_ind(pca, label="none", habillage=as.factor(all_ensg_phe$TN),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='final_pca_after_series.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_ensg_phe$subset),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='final_pca_after_studies.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_ensg_phe$datasets),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
########################################################################################
all_mir=cbind(mirExpr_tcga[common_mir,],
              mirExpr_GSE32688[common_mir,],
              mirExpr_GSE119794[common_mir,],
              mirExpr_GSE41372[common_mir,],
              mirExpr_GSE43797[common_mir,])
all_mir_phe=rbind(as.matrix(metaMatrix.MIR[,c("sample","sample_type")]),
                  as.matrix(phenoDat_GSE32688_mir[,c("newname","disease status:ch1")]),
                  as.matrix(phenoDat_GSE119794_mir[,c("newname","source_name_ch1")]),
                  as.matrix(phenoDat_GSE41372_mir[,c("newname","tissue:ch1")]),
                  as.matrix(phenoDat_GSE43797_mir[,c("newname","tissue:ch1")]))
all_mir_phe=as.data.frame(all_mir_phe)
all_mir_phe=cbind(all_mir_phe,datasets = c(rep("TCGA",nrow(metaMatrix.MIR)),
                                           rep("GSE32688",nrow(phenoDat_GSE32688_mir)),
                                           rep("GSE119794",nrow(phenoDat_GSE119794_mir)),
                                           rep("GSE41372",nrow(phenoDat_GSE41372_mir)),
                                           rep("GSE43797",nrow(phenoDat_GSE43797_mir))
))
all_mir_phe$TN=ifelse(all_mir_phe$sample_type %in% c("PrimaryTumor","pancreatic cancer","PC","pancreatic ductal adenocarcinoma","ductal adenocarcinoma (PCA)"),"Tumor",ifelse(all_mir_phe$sample_type %in% c("SolidTissueNormal","non-malignant pancreas","Normal","normal pancreas","non-neoplastic pancreas"),"Normal","wrong"))
all_mir_phe$subset=paste0(all_mir_phe$datasets,"_",all_mir_phe$TN)
all_mir=all_mir[,all_mir_phe$sample]                      

table(all_mir_phe$sample_type)
#as.matrix可以rbind 名字同原来已有的名字不相对

##示例
if(!require(factoextra))install.packages("factoextra")
if(!require(FactoMineR))install.packages("FactoMineR")
library(factoextra)
library(FactoMineR)


probesetvar = apply(all_mir, 1, var)  #表达变化大的基因
ord = order(probesetvar, decreasing=TRUE)[1:50]   #前200个基因，或者更多
pca = prcomp(t(all_mir[ord,]), scale=TRUE)
ss=summary(pca)
# 样本聚类
pdf(file='pca_before_TN.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_mir_phe$TN),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='pca_before_series.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_mir_phe$subset),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='pca_before_studies.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_mir_phe$datasets),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
# 4sigma 0.999936657516334
# 6sigma 0.999999998026825

# ,mod = NULL,这里的含义是协变量，告诉其本身有内在的分组，如果不删除tcga里癌旁正常，可以留着。 一定是先校正后删除样本 20201020
library(sva)
mod = model.matrix(~as.factor(TN), data=all_mir_phe)
all_mir=ComBat(dat=all_mir,batch=all_mir_phe$datasets,mod=mod)
# rnaExpr=as.data.frame(rnaExpr_quant) 

probesetvar = apply(all_mir, 1, var)  #表达变化大的基因
ord = order(probesetvar, decreasing=TRUE)[1:50]   #前200个基因，或者更多
pca = prcomp(t(all_mir[ord,]), scale=TRUE)
ss=summary(pca)
pdf(file='pca_after_TN.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_mir_phe$TN),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='pca_after_series.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_mir_phe$subset),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='pca_after_studies.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_mir_phe$datasets),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
#Find outlier
metaMatrix.MIR=subset(metaMatrix.MIR,!sample %in%c("TCGA-F2-6880-01","TCGA-H6-A45N-11","TCGA-2J-AABP-01"))
mirExpr_tcga=mirExpr_tcga[,metaMatrix.MIR$sample]
phenoDat_GSE32688_mir=subset(phenoDat_GSE32688_mir,!newname %in%c("humanPDAC24"))
mirExpr_GSE32688=mirExpr_GSE32688[,phenoDat_GSE32688_mir$newname]                     
#
all_mir=cbind(mirExpr_tcga[common_mir,],
              mirExpr_GSE32688[common_mir,],
              mirExpr_GSE119794[common_mir,],
              mirExpr_GSE41372[common_mir,],
              mirExpr_GSE43797[common_mir,])

all_mir_phe=rbind(as.matrix(metaMatrix.MIR[,c("sample","sample_type")]),
                  as.matrix(phenoDat_GSE32688_mir[,c("newname","disease status:ch1")]),
                  as.matrix(phenoDat_GSE119794_mir[,c("newname","source_name_ch1")]),
                  as.matrix(phenoDat_GSE41372_mir[,c("newname","tissue:ch1")]),
                  as.matrix(phenoDat_GSE43797_mir[,c("newname","tissue:ch1")]))
all_mir_phe=as.data.frame(all_mir_phe)
all_mir_phe=cbind(all_mir_phe,datasets = c(rep("TCGA",nrow(metaMatrix.MIR)),
                                           rep("GSE32688",nrow(phenoDat_GSE32688_mir)),
                                           rep("GSE119794",nrow(phenoDat_GSE119794_mir)),
                                           rep("GSE41372",nrow(phenoDat_GSE41372_mir)),
                                           rep("GSE43797",nrow(phenoDat_GSE43797_mir))
))


all_mir_phe$TN=ifelse(all_mir_phe$sample_type %in% c("PrimaryTumor","pancreatic cancer","PC","pancreatic ductal adenocarcinoma","ductal adenocarcinoma (PCA)"),"Tumor",ifelse(all_mir_phe$sample_type %in% c("SolidTissueNormal","non-malignant pancreas","Normal","normal pancreas","non-neoplastic pancreas"),"Normal","wrong"))
all_mir_phe$subset=paste0(all_mir_phe$datasets,"_",all_mir_phe$TN)
all_mir=all_mir[,all_mir_phe$sample]                      

# ,mod = NULL,这里的含义是协变量，告诉其本身有内在的分组，如果不删除tcga里癌旁正常，可以留着。 一定是先校正后删除样本 20201020
library(sva)
mod = model.matrix(~as.factor(TN), data=all_mir_phe)
all_mir=ComBat(dat=all_mir,batch=all_mir_phe$datasets,mod=mod)
# rnaExpr=as.data.frame(rnaExpr_quant) 

probesetvar = apply(all_mir, 1, var)  #表达变化大的基因
ord = order(probesetvar, decreasing=TRUE)[1:50]   #前200个基因，或者更多
pca = prcomp(t(all_mir[ord,]), scale=TRUE)
ss=summary(pca)
pdf(file='final_pca_after_TN.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_mir_phe$TN),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='final_pca_after_series.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_mir_phe$subset),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()
pdf(file='final_pca_after_studies.pdf',height = 8,width = 8)
fviz_pca_ind(pca, label="none", habillage=as.factor(all_mir_phe$datasets),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.999936657516334, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25)) #+ coord_fixed(1)
dev.off()

# rnaExpr_tcga=all_ensg
# rnaExpr_tcga,metaMatrix.RNA
# rnaExpr_GSE32688,phenoDat_GSE32688_rna
# rnaExpr_GSE119794,phenoDat_GSE119794_rna
# rnaExpr_GSE41372,phenoDat_GSE41372_rna
# rnaExpr_GSE43797,phenoDat_GSE43797_rna
# 
# all_ensg_phe
# all_mir
# all_mir_phe
selectTrait <- function(datasets,phenoDat,coln,traitV,trait){
  new_phenoDat=subset(phenoDat,phenoDat[[traitV]]==trait)
  datasets=as.data.frame(datasets)
  new_datasets=subset(datasets,select=colnames(datasets) %in% new_phenoDat[[coln]])
  return(list(new_phenoDat,new_datasets))
}
mrnaALL=list(selectTrait(all_ensg,all_ensg_phe,"sample","datasets","TCGA"),
             selectTrait(all_ensg,all_ensg_phe,"sample","datasets","GSE32688"),
             selectTrait(all_ensg,all_ensg_phe,"sample","datasets","GSE119794"),
             selectTrait(all_ensg,all_ensg_phe,"sample","datasets","GSE41372"),
             selectTrait(all_ensg,all_ensg_phe,"sample","datasets","GSE43797"))
mirnaALL=list(selectTrait(all_mir,all_mir_phe,"sample","datasets","TCGA"),
              selectTrait(all_mir,all_mir_phe,"sample","datasets","GSE32688"),
              selectTrait(all_mir,all_mir_phe,"sample","datasets","GSE119794"),
              selectTrait(all_mir,all_mir_phe,"sample","datasets","GSE41372"),
              selectTrait(all_mir,all_mir_phe,"sample","datasets","GSE43797"))


all_ensg_pheT=subset(all_ensg_phe,TN %in% "Tumor")
all_mir_pheT=subset(all_mir_phe,TN %in% "Tumor")
mrnaT=list(selectTrait(all_ensg,all_ensg_pheT,"sample","datasets","TCGA"),
           selectTrait(all_ensg,all_ensg_pheT,"sample","datasets","GSE32688"),
           selectTrait(all_ensg,all_ensg_pheT,"sample","datasets","GSE119794"),
           selectTrait(all_ensg,all_ensg_pheT,"sample","datasets","GSE41372"),
           selectTrait(all_ensg,all_ensg_pheT,"sample","datasets","GSE43797"))
mirnaT=list(selectTrait(all_mir,all_mir_pheT,"sample","datasets","TCGA"),
            selectTrait(all_mir,all_mir_pheT,"sample","datasets","GSE32688"),
            selectTrait(all_mir,all_mir_pheT,"sample","datasets","GSE119794"),
            selectTrait(all_mir,all_mir_pheT,"sample","datasets","GSE41372"),
            selectTrait(all_mir,all_mir_pheT,"sample","datasets","GSE43797"))
all_ensg_pheN=subset(all_ensg_phe,TN %in% "Normal")
all_mir_pheN=subset(all_mir_phe,TN %in% "Normal")
mrnaN=list(selectTrait(all_ensg,all_ensg_pheN,"sample","datasets","TCGA"),
           selectTrait(all_ensg,all_ensg_pheN,"sample","datasets","GSE32688"),
           selectTrait(all_ensg,all_ensg_pheN,"sample","datasets","GSE119794"),
           selectTrait(all_ensg,all_ensg_pheN,"sample","datasets","GSE41372"),
           selectTrait(all_ensg,all_ensg_pheN,"sample","datasets","GSE43797"))
mirnaN=list(selectTrait(all_mir,all_mir_pheN,"sample","datasets","TCGA"),
            selectTrait(all_mir,all_mir_pheN,"sample","datasets","GSE32688"),
            selectTrait(all_mir,all_mir_pheN,"sample","datasets","GSE119794"),
            selectTrait(all_mir,all_mir_pheN,"sample","datasets","GSE41372"),
            selectTrait(all_mir,all_mir_pheN,"sample","datasets","GSE43797"))
if(F){
  mrnaT=list(selectTrait(rnaExpr_tcga,metaMatrix.RNA,"sample","sample_type","PrimaryTumor"),
             selectTrait(rnaExpr_GSE32688,phenoDat_GSE32688_rna,"newname","disease status:ch1","pancreatic cancer"),
             selectTrait(rnaExpr_GSE119794,phenoDat_GSE119794_rna,"newname","source_name_ch1","PC"),
             selectTrait(rnaExpr_GSE41372,phenoDat_GSE41372_rna,"newname","tissue:ch1","pancreatic ductal adenocarcinoma"),
             selectTrait(rnaExpr_GSE43797,phenoDat_GSE43797_rna,"newname","tissue:ch1","ductal adenocarcinoma (PCA)"))
  
  mirnaT=list(selectTrait(mirExpr_tcga,metaMatrix.MIR,"sample","sample_type","PrimaryTumor"),
              selectTrait(mirExpr_GSE32688,phenoDat_GSE32688_mir,"newname","disease status:ch1","pancreatic cancer"),
              selectTrait(mirExpr_GSE119794,phenoDat_GSE119794_mir,"newname","source_name_ch1","PC"),
              selectTrait(mirExpr_GSE41372,phenoDat_GSE41372_mir,"newname","tissue:ch1","pancreatic ductal adenocarcinoma"),
              selectTrait(mirExpr_GSE43797,phenoDat_GSE43797_mir,"newname","tissue:ch1","ductal adenocarcinoma (PCA)"))
  
  
  mrnaN=list(selectTrait(rnaExpr_tcga,metaMatrix.RNA,"sample","sample_type","SolidTissueNormal"),
             selectTrait(rnaExpr_GSE32688,phenoDat_GSE32688_rna,"newname","disease status:ch1","non-malignant pancreas"),
             selectTrait(rnaExpr_GSE119794,phenoDat_GSE119794_rna,"newname","source_name_ch1","Normal"),
             selectTrait(rnaExpr_GSE41372,phenoDat_GSE41372_rna,"newname","tissue:ch1","normal pancreas"),
             selectTrait(rnaExpr_GSE43797,phenoDat_GSE43797_rna,"newname","tissue:ch1","non-neoplastic pancreas"))
  
  mirnaN=list(selectTrait(mirExpr_tcga,metaMatrix.MIR,"sample","sample_type","SolidTissueNormal"),
              selectTrait(mirExpr_GSE32688,phenoDat_GSE32688_mir,"newname","disease status:ch1","non-malignant pancreas"),
              selectTrait(mirExpr_GSE119794,phenoDat_GSE119794_mir,"newname","source_name_ch1","Normal"),
              selectTrait(mirExpr_GSE41372,phenoDat_GSE41372_mir,"newname","tissue:ch1","normal pancreas"),
              selectTrait(mirExpr_GSE43797,phenoDat_GSE43797_mir,"newname","tissue:ch1","non-neoplastic pancreas"))
  
  # mrnaALL=list(list(rnaExpr_tcga,metaMatrix.RNA),
  #              list(rnaExpr_GSE32688,phenoDat_GSE32688_rna),
  #              list(rnaExpr_GSE119794,phenoDat_GSE119794_rna),
  #              list(rnaExpr_GSE41372,phenoDat_GSE41372_rna),
  #              list(rnaExpr_GSE43797,phenoDat_GSE43797_rna)
  # )
  # 
  # mirnaALLL=list(list(mirExpr_tcga,metaMatrix.MIR),
  #                list(mirExpr_GSE32688,phenoDat_GSE32688_mir),
  #                list(mirExpr_GSE119794,phenoDat_GSE119794_mir),
  #                list(mirExpr_GSE41372,phenoDat_GSE41372_mir),
  #                list(mirExpr_GSE43797,phenoDat_GSE43797_mir)
  # )
}



library(doParallel) #并行处理包
cl <- makeCluster(5)#makeCluster(detectCores())
registerDoParallel(cl)

resultALL <- foreach(i=1:5) %dopar% GDCRNATools::gdcCEAnalysis(lnc        = lncrna$gene_id, 
                                                               pc          = mrna$gene_id,
                                                               deMIR       = common_mir,
                                                               lnc.targets = lnc_inter_target,#"starBase",#lnc_inter_target
                                                               pc.targets  = mrna_inter_target,#"starBase",#mrna_inter_target
                                                               rna.expr    = mrnaALL[[i]][[2]], 
                                                               mir.expr    = mirnaALL[[i]][[2]]) # apply的并行版本
resultALL1 <- foreach(i=1:5) %dopar% subset(resultALL[[i]],hyperPValue<0.05&cor>=0.0&corPValue<0.05)
resultALLlist <- foreach(i=1:5) %dopar% paste(resultALL1[[i]]$lncRNAs,resultALL1[[i]]$Genes,sep="_")

ansALL <- Reduce(intersect, resultALLlist)



reeeeeeeeeee <- foreach(i=1:5) %dopar% cbind(resultALL1[[i]],resultALLlist[[i]])
reeeeeeeeeee1 <- foreach(i=1:5) %dopar% subset(reeeeeeeeeee[[i]],reeeeeeeeeee[[i]]$`resultALLlist[[i]]` %in% ansALL)
ansALLmi=unlist(strsplit(reeeeeeeeeee1[[1]]$miRNAs,","))
ansALLmi=unique(ansALLmi)
#90个mirna
resultT <- foreach(i=1:5) %dopar% GDCRNATools::gdcCEAnalysis(lnc        = lncrna$gene_id, 
                                                             pc          = mrna$gene_id,
                                                             deMIR       = common_mir,
                                                             lnc.targets = lnc_inter_target,#"starBase",#lnc_inter_target
                                                             pc.targets  = mrna_inter_target,#"starBase",#mrna_inter_target
                                                             rna.expr    = mrnaT[[i]][[2]], 
                                                             mir.expr    = mirnaT[[i]][[2]]) # apply的并行版本
resultT1 <- foreach(i=1:5) %dopar% subset(resultT[[i]],hyperPValue<0.05&cor>=0.0&corPValue<0.05)
resultTlist <- foreach(i=1:5) %dopar% paste(resultT1[[i]]$lncRNAs,resultT1[[i]]$Genes,sep="_")
for (i in 1:length(resultNlist)){
  names(resultT1)[i] <- c("TCGA","GSE32688","GSE119794","GSE41372","GSE43797")[i]
}
write.csv(resultT1[[1]],"TCGA.csv",quote = T)
write.csv(resultT1[[2]],"GSE32688.csv",quote = T)
write.csv(resultT1[[3]],"GSE119794.csv",quote = T)
write.csv(resultT1[[4]],"GSE41372.csv",quote = T)
write.csv(resultT1[[5]],"GSE43797.csv",quote = T)

ansT <- Reduce(intersect, resultTlist)

resultN <- foreach(i=1:5) %dopar% GDCRNATools::gdcCEAnalysis(lnc        = lncrna$gene_id, 
                                                             pc          = mrna$gene_id,
                                                             deMIR       = common_mir,
                                                             lnc.targets = lnc_inter_target,#"starBase",#lnc_inter_target
                                                             pc.targets  = mrna_inter_target,#"starBase",#mrna_inter_target
                                                             rna.expr    = mrnaN[[i]][[2]], 
                                                             mir.expr    = mirnaN[[i]][[2]]) # apply的并行版本
resultN1 <- foreach(i=1:5) %dopar% subset(resultN[[i]],hyperPValue<0.05&cor>=0.0&corPValue<0.05)
resultNlist <- foreach(i=1:5) %dopar% paste(resultN1[[i]]$lncRNAs,resultN1[[i]]$Genes,sep="_")

ansN <- Reduce(intersect, resultNlist)

# resultTlist <- foreach(i=1:5) %dopar% Reduce(intersect, resultTlist)

# paste(resultT1[[i]]$lncRNAs,resultT1[[i]]$Genes,sep="_")

# as.data.frame(combn(5,3))[[10]]
stopCluster(cl)

library(gplots)
for (i in 1:length(resultNlist)){
  names(resultNlist)[i] <- c("TCGA","GSE32688","GSE119794","GSE41372","GSE43797")[i]
}
for (i in 1:length(resultTlist)){
  names(resultTlist)[i] <- c("TCGA","GSE32688","GSE119794","GSE41372","GSE43797")[i]
}
for (i in 1:length(resultALLlist)){
  names(resultALLlist)[i] <- c("TCGA","GSE32688","GSE119794","GSE41372","GSE43797")[i]
}
save(resultT,resultTlist,resultN,resultNlist,resultALL,resultALLlist,file="f2021.Rdata")
# load(file="f.Rdata")
kkkk1=gplots::venn(resultNlist)
kkkk2=gplots::venn(resultTlist)
kkkk3=gplots::venn(resultALLlist)
pdf("bo.pdf", height=10, width=9) 
par(mfrow=c(1,2))
plot(kkkk1,kkkk2,kkkk3)
dev.off() 

# outPath <- "./deg/TAD fusion/" ##输出路径
# out_fileName <- sapply(names(resultALL2),function(x){
#   paste(x, ".csv", sep='')}) ##csv格式
# out_filePath  <- sapply(out_fileName, function(x){
#   paste(outPath ,x,sep='/')}) ##输出路径名
# ##输出文件
# for(i in 1:length(resultALL2)){
#   write.csv(resultALL2[[i]], file=out_filePath[i],quote = T) 
# }
# 
# 




library(VennDiagram)
venn.diagram(resultALLlist,filename = "1Venn.tiff",height = 10000, width = 10000,resolution =1200,cex=1.5,cat.cex = 2,imagetype = "tiff",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,	cat.dist = c(0.2, 0.2, 0.15, 0.22,0.2),cat.pos = c(0,-20,180,180,30),margin = 0.02)
venn.diagram(resultTlist,filename = "2Venn.tiff",height = 10000, width = 10000,resolution =1200,cex=1.5,cat.cex = 2,imagetype = "tiff",
             fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
             alpha = 0.50,	cat.dist = c(0.2, 0.2, 0.15, 0.22,0.2),cat.pos = c(0,-20,180,180,30),margin = 0.02)
venn.diagram(resultNlist,filename = "3Venn.tiff",height = 10000, width = 10000,resolution =1200,cex=1.5,cat.cex = 2,imagetype = "tiff",
             fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
             alpha = 0.50,	cat.dist = c(0.2, 0.2, 0.15, 0.22,0.2),cat.pos = c(0,-20,180,180,30),margin = 0.02)
library(UpSetR)

library(GDCRNATools)
library(VennDiagram)#只能有五个
# resultTlist_U1=calculate.overlap(resultTlist)
resultTlist_U=get.venn.partitions(resultTlist)
resultTlist_U[,c(1:5)]=apply(resultTlist_U[,c(1:5)],2,function(x) as.numeric(x))
resultTlist_U$sum=apply(resultTlist_U[,c(1:5)],1,sum)
# 3D 热图 20211108
resultT_ce_mir=resultT[[1]]
resultT_ce_mir$miRNAs=as.character(resultT_ce_mir$miRNAs)
resultT_ce_mir=tidyr::separate_rows(resultT_ce_mir,miRNAs,sep=",")
resultT_ce_mir=unite(resultT_ce_mir,col = "e",c("lncRNAs","Genes"),remove = F)
resultTlist_U2=resultTlist_U
resultTlist_U2=tidyr::unnest_longer(resultTlist_U2,..values..)
resultTlist_frequency=merge(resultT_ce_mir,resultTlist_U2,by.x="e",by.y="..values..")
save(resultTlist_frequency,file="resultTlist_frequency.Rdata")
#
ansT3_equal <- Reduce(union,resultTlist_U[resultTlist_U$sum>2,]$..values..)
ansT4_equal <- Reduce(union,resultTlist_U[resultTlist_U$sum>3,]$..values..)
resultTlist_U=resultTlist_U[resultTlist_U$TCGA=="1",]
resultTlist_U=resultTlist_U[resultTlist_U$sum>2,]
ansT3 <- Reduce(union,resultTlist_U$..values..)
ansT4 <- Reduce(union,resultTlist_U[resultTlist_U$sum>3,]$..values..)
mT=stringr::str_split(ansT3,"_",simplify = T)[,2]
lncT=stringr::str_split(ansT3,"_",simplify = T)[,1]
resultT_ce=subset(resultT1[[1]],resultT1[[1]]$Genes%in%mT &resultT1[[1]]$lncRNAs %in%lncT)
edges <- gdcExportNetwork(ceNetwork = resultT_ce, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = resultT_ce, net = 'nodes')
write.csv(edges, file='edgesT.csv',quote = F) ### Network of Cytoscape
write.csv(nodes, file='nodesT.csv',quote = F) ### Table of Cytoscape

library(tidyverse)
resultTlist_U_csv=resultTlist_U
resultTlist_U_csv=data.frame(rep(paste0(resultTlist_U_csv$TCGA,"_",resultTlist_U_csv$GSE32688,"_",resultTlist_U_csv$GSE119794,"_",resultTlist_U_csv$GSE41372,"_",resultTlist_U_csv$GSE43797),resultTlist_U_csv$..count..),unlist(resultTlist_U_csv$..values..))
resultTlist_U_csv=cbind(resultTlist_U_csv,str_split(resultTlist_U_csv$rep.paste0.resultTlist_U_csv.TCGA..._...resultTlist_U_csv.GSE32688..,"_",simplify = T))
colnames(resultTlist_U_csv)=c("coll","ensg","TCGA","GSE32688","GSE119794","GSE41372","GSE43797")
write.csv(resultTlist_U_csv, file='resultTlist_U_csv.csv',quote = T) ### Table of Cytoscape


library(VennDiagram)#只能有五个
# resultTlist_U1=calculate.overlap(resultTlist)
resultALLlist_U=get.venn.partitions(resultALLlist)
resultALLlist_U[,c(1:5)]=apply(resultALLlist_U[,c(1:5)],2,function(x) as.numeric(x))
resultALLlist_U$sum=apply(resultALLlist_U[,c(1:5)],1,sum)
ansALL3_equal <- Reduce(union,resultALLlist_U[resultALLlist_U$sum>2,]$..values..)
ansALL4_equal <- Reduce(union,resultALLlist_U[resultALLlist_U$sum>3,]$..values..)
resultALLlist_U=resultALLlist_U[resultALLlist_U$TCGA=="1",]
resultALLlist_U=resultALLlist_U[resultALLlist_U$sum>2,]
ansALL3 <- Reduce(union,resultALLlist_U$..values..)
ansALL4 <- Reduce(union,resultALLlist_U[resultALLlist_U$sum>3,]$..values..)
mALL=stringr::str_split(ansALL3,"_",simplify = T)[,2]
lncALL=stringr::str_split(ansALL3,"_",simplify = T)[,1]
resultALL_ce=subset(resultALL1[[1]],resultALL1[[1]]$Genes%in%mALL &resultALL1[[1]]$lncRNAs %in%lncALL)
edges <- gdcExportNetwork(ceNetwork = resultALL_ce, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = resultALL_ce, net = 'nodes')
write.csv(edges, file='edgesALL.csv',quote = F) ### Network of Cytoscape
write.csv(nodes, file='nodesALL.csv',quote = F) ### Table of Cytoscape


library(VennDiagram)#只能有五个
# resultTlist_U1=calculate.overlap(resultTlist)
resultNlist_U=get.venn.partitions(resultNlist)
resultNlist_U[,c(1:5)]=apply(resultNlist_U[,c(1:5)],2,function(x) as.numeric(x))
resultNlist_U$sum=apply(resultNlist_U[,c(1:5)],1,sum)
ansN3_equal <- Reduce(union,resultNlist_U[resultNlist_U$sum>2,]$..values..)
ansN4_equal <- Reduce(union,resultNlist_U[resultNlist_U$sum>3,]$..values..)
resultNlist_U=resultNlist_U[resultNlist_U$TCGA=="1",]
resultNlist_U=resultNlist_U[resultNlist_U$sum>2,]
ansN3 <- Reduce(union,resultNlist_U$..values..)
ansN4 <- Reduce(union,resultNlist_U[resultNlist_U$sum>3,]$..values..)
mN=stringr::str_split(ansN3,"_",simplify = T)[,2]
lncN=stringr::str_split(ansN3,"_",simplify = T)[,1]
resultN_ce=subset(resultN1[[1]],resultN1[[1]]$Genes%in%mN &resultN1[[1]]$lncRNAs %in%lncN)
edges <- gdcExportNetwork(ceNetwork = resultN_ce, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = resultN_ce, net = 'nodes')
write.csv(edges, file='edgesN.csv',quote = F) ### Network of Cytoscape
write.csv(nodes, file='nodesN.csv',quote = F) ### Table of Cytoscape

save(resultALL,resultT,resultN,ansT3,ansT4,ansALL3,ansALL4, ansN3,ansN4,file="ans.Rdata")

load("C:/Users/zxh/Desktop/R/paad-tcga-gtex/kaiqiaole_v2.Rdata")
mat.out=unique(rbind(mat.out10,mat.out5,mat.out3))#mat.outaa,
mat.out$number=sapply(stringr::str_split(mat.out$X8,"\\+",simplify = F),function(x) length(x))
mat.out=subset(mat.out,number==7)
T7=unlist(strsplit(mat.out$X8,"\\+"))

resultT_ce_7=resultT_ce[resultT_ce$Genes %in%T7,]
edges <- gdcExportNetwork(ceNetwork = resultT_ce_7, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = resultT_ce_7, net = 'nodes')
write.csv(edges, file='edgesT7.csv',quote = F) ### Network of Cytoscape
write.csv(nodes, file='nodesT7.csv',quote = F) ### Table of Cytoscape

load(file="ans.Rdata")

all3=list(ansALL3_equal,ansT3_equal,ansN3_equal)
all4=list(ansALL4_equal,ansT4_equal,ansN4_equal)
for (i in 1:length(all3)){
  names(all3)[i] <- c("All","Tumor","Normal")[i]
}
for (i in 1:length(all4)){
  names(all4)[i] <- c("All","Tumor","Normal")[i]
}
library(VennDiagram)
venn.diagram(all3,filename = "all3.tiff",col = "transparent",height = 10000, width = 10000,resolution =1200, imagetype = "tiff",  cex = 2,  cat.default.pos = "text",  cat.col = c("darkred", "darkblue", "darkgreen"),  cat.cex = 2,  cat.dist = c(0.28, 0.06, 0.03),fill = c("dodgerblue", "goldenrod1", "darkorange1"),alpha = 0.50,	cat.pos = 0,margin = 0.02)
venn.diagram(all4,filename = "all4.tiff",col = "transparent",height = 10000, width = 10000,resolution =1200, imagetype = "tiff",  cex = 2,  cat.default.pos = "text",  cat.col = c("darkred", "darkblue", "darkgreen"),  cat.cex = 2,  cat.dist = c(0.06, 0.06, 0.03),fill = c("dodgerblue", "goldenrod1", "darkorange1"),alpha = 0.50,	cat.pos = 0,margin = 0.02)


pdf(file='upset_mrna_1.pdf',height = 8,width = 8,onefile = T)

# for (i in 1:length(resultTlist)){
#   names(resultTlist)[i] <- paste("data_", i, sep = "")
# }

# names(resultTlist2)<-c("n",names(resultTlist[-1]))
upset(fromList(resultTlist), #输入的数据集
      # sets = c("a","b"),#,"Survival"
      sets = c("TCGA","GSE32688","GSE119794","GSE41372","GSE43797"),
      # sets = c("PITA", "miRanda","TargetScan"),
      # nsets = 5, #想要可视化的数据集数量,也可以用sets选项自定义
      nintersects = 12, #要绘制的交点数量 2^n
      keep.order = T, #是否保持输入顺序,否则按各集合大小排序
      matrix.color='forestgreen',
      main.bar.color = 'forestgreen', #主条图的颜色
      mainbar.y.label = 'Number', #y轴标签
      sets.bar.color = 'red', #设置集合条图的颜色
      sets.x.label = 'Number of elements', #设置集合调图的坐标轴标签
      mb.ratio = c(0.7,0.3), #条形图点图的比例
      order.by = c("freq","degree"), #交集如何排序,
      decreasing = c(TRUE,TRUE),
      text.scale =1.9,
      # 以上排序是否降序,FALSE
      # boxplot.summary = c('ReleaseDate','Comedy') #添加箱线图展示数据分布,最多展示两个数据集
      queries = list(list(query = intersects,
                          params = list("TCGA","GSE32688","GSE119794","GSE41372","GSE43797"),#,"Survival"
                          color = 'black',
                          active = T)
                     # ,
      #                list(query = intersects,
      #                     params = list("Gene_down","Gene_predict", "Gene_WGCNA"),#,"Survival"
      #                     color = 'green',
      #                     active = T)
      )
)

dev.off()
# myfunc <- function(row, release, rating){  
#   newdata <- (row["ReleaseDate"]%in%release)&(row["AvgRating"]>rating)  
# }  
# query=myfunc
#取出高中文表达的mirna，而不是找差异的miRNA，至少这篇不是#b
miRNA_1=mirCounts[order(apply(mirCounts,1,mad), decreasing = T)[500:0],]
mirExpr<- gdcVoomNormalization_modify(counts = mirCounts, filter = F)#limma filter
rnaExpr<- gdcVoomNormalization_modify(counts = rnaCounts, filter = F)#limma filter

ceOutput <- gdcCEAnalysis(lnc         = lncrna$gene_id, 
                          pc          = mrna$gene_id,
                          deMIR       = rownames(miRNA_1),
                          lnc.targets = lnc_inter_target,#"starBase",#lnc_inter_target
                          pc.targets  = mrna_inter_target,#"starBase",#mrna_inter_target
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)#超几何分布及p值

mirExpr <- gdcVoomNormalization_modify(counts = probes_expr_GSE119794_mir)#limma filter
rnaExpr<- gdcVoomNormalization_modify(counts = probes_expr_GSE119794_rna)#limma filter

# 是不是要过滤一下，不同数据集分别过滤吧，但不需要肿瘤和正常要分开过滤(似乎不需要)
#或者直接来干干？？
DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma',#edgeR
                        filter =T)#默认是true，不筛选则选择false
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all', fc = 0, pval = 1)#
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 1)#
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 1)#其实放到0.05也可以
LNC_0_F <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 100)#
PC_0_F <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 100)

#
probes_expr_GSE119794_mir=as.data.frame(probes_expr_GSE119794_mir)
probes_expr_GSE119794_rna=as.data.frame(probes_expr_GSE119794_rna)
forthetree=updateName(probes_expr_GSE119794_rna)
table(forthetree$gene_biotype)
probes_expr_GSE119794_rna$ALIAS=rownames(probes_expr_GSE119794_rna)
probes_expr_GSE119794_rna=merge(probes_expr_GSE119794_rna,forthetree,by="ALIAS")
rownames(probes_expr_GSE119794_rna)=probes_expr_GSE119794_rna$ENSEMBL
probes_expr_GSE119794_rna=probes_expr_GSE119794_rna[,-which(colnames(probes_expr_GSE119794_rna) %in% c("ALIAS","ENSEMBL","gene_biotype"))]

probes_expr_GSE119794_rna=unique(probes_expr_GSE119794_rna)
unique(probes_expr_GSE119794_rna$SYMBOL)

# forthetree3=forthetree[forthetree$gene_biotype %in% c("miRNA"),]
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC),
                          deMIR       = rownames(miRNA_1),
                          lnc.targets = lnc_inter_target,
                          pc.targets  = mrna_inter_target,
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)#超几何分布及p值

mirExpr <- gdcVoomNormalization_modify(counts = probes_expr_GSE119794_mir)#limma filter
rnaExpr<- gdcVoomNormalization_modify(counts = probes_expr_GSE119794_rna)#limma filter
DEGAll <- gdcDEAnalysis(counts     = probes_expr_GSE119794_rna,
                        group      = TcgaTargetGTEX_phenotype_Pancreas624$sample_type,
                        comparison = 'PrimaryTumor-SolidTissueNormal',
                        method     = 'limma',#edgeR
                        filter =T)#默认是true，不筛选则选择false
# deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all', fc = 1.000, pval = 0.05)#
# deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 1.000, pval = 0.05)#
# dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 1.000, pval = 0.05)#其实放到0.05也可以
# LNC_0_F <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 100)#
# PC_0_F <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 100)
rnaExpr_GSE41372["ENSG00000112576",]


write.csv(metaMatrix.RNA,file="fff.csv",quote =T) 