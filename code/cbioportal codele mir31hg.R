rm(list=ls())
options(stringsAsFactors = F)
gc()
library(GDCRNATools)

CDKN2Atxt=read.table(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/hic scripts/ccle 4 gene/cdkn2a.txt",sep = "\t",header = T)
CDKN2Btxt=read.csv(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/hic scripts/ccle 4 gene/cdkn2b.txt",sep = "\t",header = T)
MTAPtxt=read.csv(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/hic scripts/ccle 4 gene/mtap.txt",sep = "\t",header = T)
colnames(CDKN2Atxt)
codelete=merge(CDKN2Atxt,CDKN2Btxt,by="Sample.Id")
colnames(codelete)
codelete=merge(codelete,MTAPtxt,by="Sample.Id")
codelete=tidyr::unite(codelete,"cnv",CDKN2A..Copy.Number.Alterations,CDKN2B..Copy.Number.Alterations,MTAP..Copy.Number.Alterations,sep=" ",remove = F)
table(codelete$cnv)
codelete=subset(codelete,cnv %in% c("Deep Deletion Deep Deletion Deep Deletion","Diploid Diploid Diploid"),select = c("Sample.Id","cnv"))
MIR31HGtxt=read.csv(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/hic scripts/ccle 4 gene/mir31hg.txt",sep = "\t",header = T)
codelete=merge(codelete,MIR31HGtxt,by="Sample.Id")
codelete$cnv=ifelse(codelete$cnv=="Diploid Diploid Diploid","Normal","Co-DEL")
table(codelete$Copy.Number.Alterations)
table(codelete$MIR31HG..Copy.Number.Alterations)
codelete=subset(codelete,!(cnv %in% "Normal"&MIR31HG..Copy.Number.Alterations != "Diploid"))
codelete$cnv4=ifelse(codelete$cnv=="Normal","No changes",codelete$MIR31HG..Copy.Number.Alterations)
box=codelete
box$number1=as.numeric(box$MIR31HG..mRNA.expression..RNA.Seq.RPKM...log2.value...1.. )
box=as.data.frame(box)

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
tapply(box$MIR31HG..mRNA.expression..RNA.Seq.RPKM...log2.value...1..,box$cnv4,shapiro.test)#不正态
bartlett.test(ENSG00000171889~cnv,data=box)#方差齐性
pdf('Figure A.Boxplot of MIR31HG baed on status of codelete ccle.pdf',width = 8,height=4)
if(require('ggpubr')){
  library(ggpubr)
  # google search : ggpubr boxplot add p-value
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  p <- ggboxplot(box, x = "cnv4", y = "number1",group="cnv4",
                 color = "cnv4", palette=c(1,2,3,4),order=c("No changes","Deep Deletion","Diploid","Amplification"  ), font.label = list(size = 1, color = "black"),outlier.shape = 19#"#00AFBB","#E7B800",
  )+#add = "jitter"
    # scale_color_manual(values=c(1,2,3,4))  +#"#E69F00", "#56B4E9"
    coord_cartesian(ylim = c(0.0000, 30))+
    # theme(legend.position="none")+
    labs(title="",x="MIR31HG", y = "Expression")#+#tag="ENSG00000000003"

    # geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.2,fill = c("#999999"))
  #  Add p-value
  p + stat_compare_means(aes(x = cnv4, y = number1),method = "kruskal.test", label.y = 28)+ #备注3
    stat_compare_means(comparisons = list(c("No changes","Deep Deletion"),c("No changes","Diploid"),c("No changes","Amplification"),c("Amplification","Diploid")),method = "wilcox.test",label.y = c(5,7,12,9))+ #备注3
    stat_summary(aes(x = cnv4, y = number1),fun= "mean", geom = "point",position = position_dodge(0.75), shape = 23, size = 1, fill = "pink")
  
}
dev.off()
write.csv(box,file="co-CDKN2A-CDKN2B-MTAP-ccle.csv",quote=T)
