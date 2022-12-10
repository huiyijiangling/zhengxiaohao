rm(list=ls())
options(stringsAsFactors = F)
gc()
library(GDCRNATools)
#最重要的ajust是TSS67,center2627
load("TcgaTargetGTEX_phenotype_Pancreas.Rdata")
load("PAAD_GDCRNATOOLS.Rdata")
load("rnaCountslog21_Pancreas.Rdata")
load("rnatpmlog2001_Pancreas.Rdata")
load("rnainter_homo_sig_cerna_all.Rdata")
# load("metaMatrix.RNA.ESCA.Rdata") #没有
#TCGA-SW-A7EB-01 有信息但是没有count数据
# metaMatrix.RNA.PAAD=metaMatrix.RNA
# metaMatrix.RNA=rbind(metaMatrix.RNA.PAAD,metaMatrix.RNA.ESCA)#写错了
#数据处理RSEM expected count可以当做raw count处理，并且被TMM+voom，
#RSEM 但是找DE时，ENseq>DEseq2>≈edgeR
#z注意删除 注意质控 在运行完之后回头删除
table(TcgaTargetGTEX_phenotype_Pancreas$sample_type)
rownames(TcgaTargetGTEX_phenotype_Pancreas[TcgaTargetGTEX_phenotype_Pancreas$sample_type=="WRONG",])#"TCGA-HZ-A9TJ-06"
TcgaTargetGTEX_phenotype_Pancreas624=TcgaTargetGTEX_phenotype_Pancreas[-which(rownames(TcgaTargetGTEX_phenotype_Pancreas) %in% c("TCGA-HZ-A9TJ-06")),]
table(TcgaTargetGTEX_phenotype_Pancreas624$sample_type)
rownames(rnatpmlog2001_Pancreas)=substr(rownames(rnatpmlog2001_Pancreas),1,15)
rownames(rnaCountslog21_Pancreas)=substr(rownames(rnaCountslog21_Pancreas),1,15)

rnatpmlog2001_Pancreas=subset(rnatpmlog2001_Pancreas,select=rownames(TcgaTargetGTEX_phenotype_Pancreas624))
rnaCountslog21_Pancreas=subset(rnaCountslog21_Pancreas,select=rownames(TcgaTargetGTEX_phenotype_Pancreas624))
rnatpm_Pancreas=2^rnatpmlog2001_Pancreas-0.001#625
rnaCounts_Pancreas=2^rnaCountslog21_Pancreas-1#624
TcgaTargetGTEX_phenotype_Pancreas624codelete1=subset(TcgaTargetGTEX_phenotype_Pancreas624,TcgaTargetGTEX_phenotype_Pancreas624$sample_type %in% "SolidTissueNormal")
TcgaTargetGTEX_phenotype_Pancreas624codelete1=subset(TcgaTargetGTEX_phenotype_Pancreas624codelete1,TcgaTargetGTEX_phenotype_Pancreas624codelete1$`_study` %in% "TCGA")
TcgaTargetGTEX_phenotype_Pancreas624codelete2=subset(TcgaTargetGTEX_phenotype_Pancreas624,TcgaTargetGTEX_phenotype_Pancreas624$sample_type %in% "PrimaryTumor")

CDKN2Atxt=read.csv(file="CDKN2A.csv",skip = 0)
CDKN2Btxt=read.csv(file="CDKN2B.csv",skip = 0)
MTAPtxt=read.csv(file="MTAP.csv",skip = 0)
colnames(CDKN2Atxt)
codelete=merge(CDKN2Atxt,CDKN2Btxt,by="Sample.Id")
colnames(codelete)
codelete=merge(codelete,MTAPtxt,by="Sample.Id")
codelete=subset(codelete,MTAP..Putative.copy.number.alterations.from.GISTIC %in% "Deep Deletion")
codelete=subset(codelete,CDKN2A..Putative.copy.number.alterations.from.GISTIC %in% "Deep Deletion")
codelete=subset(codelete,CDKN2B..Putative.copy.number.alterations.from.GISTIC %in% "Deep Deletion")
codelete$cnv="Deep Deletion"
TcgaTargetGTEX_phenotype_Pancreas624codelete2=subset(TcgaTargetGTEX_phenotype_Pancreas624codelete2,rownames(TcgaTargetGTEX_phenotype_Pancreas624codelete2) %in%
                                                       codelete$Sample.Id)
TcgaTargetGTEX_phenotype_Pancreas624codelete2$Sample.Id=rownames(TcgaTargetGTEX_phenotype_Pancreas624codelete2)
TcgaTargetGTEX_phenotype_Pancreas624codelete2=merge(TcgaTargetGTEX_phenotype_Pancreas624codelete2,codelete)
rownames(TcgaTargetGTEX_phenotype_Pancreas624codelete2)=TcgaTargetGTEX_phenotype_Pancreas624codelete2$Sample.Id
TcgaTargetGTEX_phenotype_Pancreas624codelete2$cnv=TcgaTargetGTEX_phenotype_Pancreas624codelete2$cnv
TcgaTargetGTEX_phenotype_Pancreas624codelete2=TcgaTargetGTEX_phenotype_Pancreas624codelete2[,-1]
TcgaTargetGTEX_phenotype_Pancreas624codelete2=TcgaTargetGTEX_phenotype_Pancreas624codelete2[,-(8:19)]
TcgaTargetGTEX_phenotype_Pancreas624codelete1$cnv="Normal"
TcgaTargetGTEX_phenotype_Pancreas624=rbind(TcgaTargetGTEX_phenotype_Pancreas624codelete1,TcgaTargetGTEX_phenotype_Pancreas624codelete2)

rnaCountslog21_Pancreas=subset(rnaCountslog21_Pancreas,select = rownames(TcgaTargetGTEX_phenotype_Pancreas624))                                                 
rnaCounts_Pancreas=subset(rnaCounts_Pancreas,select = rownames(TcgaTargetGTEX_phenotype_Pancreas624))
rnaCounts_Pancreas=round(rnaCounts_Pancreas,0)
rnatpmlog2001_Pancreas=subset(rnatpmlog2001_Pancreas,select = rownames(TcgaTargetGTEX_phenotype_Pancreas624))     
rnatpm_Pancreas=subset(rnatpm_Pancreas,select = rownames(TcgaTargetGTEX_phenotype_Pancreas624))     
##################################################tcga
box=TcgaTargetGTEX_phenotype_Pancreas624
# box$cnv=ifelse(TcgaTargetGTEX_phenotype_Pancreas624$cnv=="PrimaryTumor",'tumor','normal')
box=cbind(box,t(rnaCounts_Pancreas["ENSG00000171889",]))
box=box[,-(1:7)]
write.csv(box,file="co-CDKN2A-CDKN2B-MTAP-count4para.csv",quote=T)
box=TcgaTargetGTEX_phenotype_Pancreas624
# box$cnv=ifelse(TcgaTargetGTEX_phenotype_Pancreas624$cnv=="PrimaryTumor",'tumor','normal')
box=cbind(box,t(rnaCountslog21_Pancreas["ENSG00000171889",]))
box=box[,-(1:7)]
write.csv(box,file="co-CDKN2A-CDKN2B-MTAP-log21-count4para.csv",quote=T)
#
box=TcgaTargetGTEX_phenotype_Pancreas624
# box$cnv=ifelse(TcgaTargetGTEX_phenotype_Pancreas624$cnv=="PrimaryTumor",'tumor','normal')
box=cbind(box,t(rnatpmlog2001_Pancreas["ENSG00000171889",]))
box=box[,-(1:7)]
write.csv(box,file="co-CDKN2A-CDKN2B-MTAP-log2001-tpm4para.csv",quote=T)

box=TcgaTargetGTEX_phenotype_Pancreas624
# box$cnv=ifelse(TcgaTargetGTEX_phenotype_Pancreas624$cnv=="PrimaryTumor",'tumor','normal')
box=cbind(box,t(rnatpm_Pancreas["ENSG00000171889",]))
box=box[,-(1:7)]
write.csv(box,file="co-CDKN2A-CDKN2B-MTAP-tpm4para.csv",quote=T)
# pdf('Figure A.Boxplot of codelete 4para.pdf',width = 8,height=4)
# if(require('ggpubr')){
#   library(ggpubr)
#   # google search : ggpubr boxplot add p-value
#   # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
#   p <- ggboxplot(box, x = "cnv", y = "ENSG00000147883",
#                  color = "cnv", palette=c(1,2,3,4,5),order = c("Normal","Deep Deletion","Shallow Deletion" ,"Diploid","Gain"),  font.label = list(size = 1, color = "black"),#"#00AFBB","#E7B800"
#   )+#add = "jitter"
#     scale_color_manual(values=c(1,2,3,4,5))  +#"#E69F00", "#56B4E9"
#     theme(legend.position="none")+
#     labs(title="",x="codelete", y = "Expression")+#tag="ENSG00000000003"
#     geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.2,fill = c("#999999"))
#   #  Add p-value
#   p + stat_compare_means(label.y.npc="bottom",label = "p.format")
# }
# dev.off()

# pdf('Figure A.Boxplot of MLLT3 baed on status of codelete 4para.pdf',width = 8,height=4)
# if(require('ggpubr')){
#   library(ggpubr)
#   # google search : ggpubr boxplot add p-value
#   # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
#   p <- ggboxplot(box, x = "cnv", y = "ENSG00000171843",
#                  color = "cnv", palette=c(1,2,3,4,5),order = c("Normal","Deep Deletion","Shallow Deletion" ,"Diploid","Gain"),  font.label = list(size = 1, color = "black"),#"#00AFBB","#E7B800"
#   )+#add = "jitter"
#     scale_color_manual(values=c(1,2,3,4,5))  +#"#E69F00", "#56B4E9"
#     
#     theme(legend.position="none")+
#     labs(title="",x="MLLT3", y = "Expression")+#tag="ENSG00000000003"
#     geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.2,fill = c("#999999"))
#   #  Add p-value
#   p + stat_compare_means(label.y.npc="bottom",label = "p.format")
# }
# dev.off()
# ks.test(BxPC3$len)
shapiro.test#3-5000
lillie.test#>5000
library(nortest)
tapply(box$ENSG00000171889,box$cnv,shapiro.test)#不正态
bartlett.test(ENSG00000171889~cnv,data=box)#方差齐性
pdf('Figure A.Boxplot of MIR31HG baed on status of codelete 4para.pdf',width = 8,height=4)
if(require('ggpubr')){
  library(ggpubr)
  # google search : ggpubr boxplot add p-value
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  p <- ggboxplot(box, x = "cnv", y = "ENSG00000171889",
                 color = "cnv", palette=c(1,2,3,4,5),order = c("Normal","Deep Deletion","Shallow Deletion" ,"Diploid","Gain"),  font.label = list(size = 1, color = "black"),#"#00AFBB","#E7B800"
  )+#add = "jitter"
    scale_color_manual(values=c(1,2,3,4,5))  +#"#E69F00", "#56B4E9"
    
    theme(legend.position="none")+
    labs(title="",x="MIR31HG", y = "Expression")+#tag="ENSG00000000003"
    geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.2,fill = c("#999999"))
  #  Add p-value
  p + stat_compare_means(label.y.npc="top",method="wilcox.test")#,label = "p.format"
}
dev.off()
?stat_compare_means
