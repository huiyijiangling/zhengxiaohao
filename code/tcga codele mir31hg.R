rm(list=ls())
options(stringsAsFactors = F)
gc()
library(GDCRNATools)

CDKN2Atxt=read.table(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/hic scripts/tcga 4 gene/cdkn2a.txt",sep = "\t",header = T)
CDKN2Btxt=read.csv(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/hic scripts/tcga 4 gene/cdkn2b.txt",sep = "\t",header = T)
MTAPtxt=read.csv(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/hic scripts/tcga 4 gene/mtap.txt",sep = "\t",header = T)
colnames(CDKN2Atxt)
codelete=merge(CDKN2Atxt,CDKN2Btxt,by="Sample.Id")
colnames(codelete)
codelete=merge(codelete,MTAPtxt,by="Sample.Id")
codelete=tidyr::unite(codelete,"cnv",CDKN2A..Putative.copy.number.alterations.from.GISTIC,CDKN2B..Putative.copy.number.alterations.from.GISTIC,MTAP..Putative.copy.number.alterations.from.GISTIC,sep=" ",remove = F)
table(codelete$cnv)
codelete=subset(codelete,cnv %in% c("Deep Deletion Deep Deletion Deep Deletion","Diploid Diploid Diploid","Shallow Deletion Shallow Deletion Shallow Deletion"),select = c("Sample.Id","cnv"))
MIR31HGtxt=read.csv(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/hic scripts/tcga 4 gene/mir31hg.txt",sep = "\t",header = T)
codelete=merge(codelete,MIR31HGtxt,by="Sample.Id")
codelete$cnv=ifelse(codelete$cnv=="Diploid Diploid Diploid","No changes",
                    ifelse(codelete$cnv=="Shallow Deletion Shallow Deletion Shallow Deletion","Shallow Deletion","Deep Deletion"))
table(codelete$Copy.Number.Alterations)
table(codelete$MIR31HG..Putative.copy.number.alterations.from.GISTIC)
# codelete=subset(codelete,!(cnv %in% "No changes"&MIR31HG..Putative.copy.number.alterations.from.GISTIC != "Diploid"))#9
codelete$cnv4=ifelse(codelete$cnv!="Deep Deletion",codelete$cnv,codelete$Copy.Number.Alterations)
box=codelete
box$number1=as.numeric(box$MIR31HG..mRNA.Expression..RSEM..Batch.normalized.from.Illumina.HiSeq_RNASeqV2...log2.value...1..)
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
# table(box$cnv4)
tapply(box$MIR31HG..mRNA.Expression..RSEM..Batch.normalized.from.Illumina.HiSeq_RNASeqV2...log2.value...1..,box$cnv4,lillie.test)#不正态
bartlett.test(ENSG00000171889~cnv,data=box)#方差齐性
pdf('Figure A.Boxplot of MIR31HG baed on status of codelete ccle.pdf',width = 8,height=8)
if(require('ggpubr')){
  library(ggpubr)
  # google search : ggpubr boxplot add p-value
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  p <- ggboxplot(box, x = "cnv4", y = "number1",group="cnv4",
                 color = "cnv4",  font.label = list(size = 1, color = "black"),palette=c(1,2,3,4,5,6,7),order=c("No changes","Shallow Deletion","MIR31HG: Deep Deletion","MIR31HG: Shallow Deletion","MIR31HG: Diploid","MIR31HG: Gain","MIR31HG: Amplification"),outlier.shape = NA#"#00AFBB","#E7B800",
  )+#add = "jitter"
    # scale_color_manual(values=c(1,2,3,4))  +#"#E69F00", "#56B4E9"
    coord_cartesian(ylim = c(0.0000, 500))+
    # theme(legend.position="none")+
    labs(title="",x="MIR31HG", y = "Expression")#+#tag="ENSG00000000003"
  
  # geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.2,fill = c("#999999"))
  #  Add p-value
  p + stat_compare_means(aes(x = cnv4, y = number1),method = "kruskal.test", label.y = 500)+ #备注3
    stat_compare_means(comparisons = list(c("No changes","Shallow Deletion"),c("No changes","MIR31HG: Deep Deletion"),c("No changes","MIR31HG: Shallow Deletion"),c("No changes","MIR31HG: Diploid"),c("No changes","MIR31HG: Gain"),c("No changes","MIR31HG: Amplification")),method = "wilcox.test",label.y = as.numeric(seq(100,850,by=20)))+ #备注3
    stat_summary(aes(x = cnv4, y = number1),fun= "mean", geom = "point",position = position_dodge(0.75), shape = 23, size = 1, fill = "pink")

}
dev.off()
testpn
write.csv(box,file="co-CDKN2A-CDKN2B-MTAP-tcga.csv",quote=T)

cellall=readxl::read_xlsx("C:/Users/zxh/Desktop/R/original cellPAAD.xlsx",sheet=3)
table(duplicated(substr(box$Sample.Id,1,12)))#false
box$bcr_patient_barcode=substr(box$Sample.Id,1,12)
box_surv=merge(box,cellall,by.x="bcr_patient_barcode",by.y = "bcr_patient_barcode")
dat=box_surv


dat$event=as.numeric(dat$OS)#dat$event=as.numeric(dat$DFS)#dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$OS.time)#dat$time=as.numeric(dat$DFS.time)#dat$time=as.numeric(dat$time)
#
fit <- survfit(Surv(time, event) ~cnv4, data = dat,type  = "kaplan-meier",conf.type = "log-log")
pdf("os.pdf",height = 8,width = 16,onefile = F)
ggsurvplot(fit, data = dat, risk.table = TRUE, conf.int = F, pval=TRUE, break.time.by = 365)
dev.off()
#
dat$event=as.numeric(dat$DFI)#dat$event=as.numeric(dat$DFS)#dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$DFI.time)#dat$time=as.numeric(dat$DFS.time)#dat$time=as.numeric(dat$time)
fit <- survfit(Surv(time, event) ~cnv4, data = dat,type  = "kaplan-meier",conf.type = "log-log")
pdf("dfi.pdf",height = 8,width = 16,onefile = F)
ggsurvplot(fit, data = dat, risk.table = TRUE, conf.int = F, pval=TRUE, break.time.by = 365)
dev.off()
#
dat$event=as.numeric(dat$PFI)#dat$event=as.numeric(dat$DFS)#dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$PFI.time)#dat$time=as.numeric(dat$DFS.time)#dat$time=as.numeric(dat$time)
fit <- survfit(Surv(time, event) ~cnv4, data = dat,type  = "kaplan-meier",conf.type = "log-log")
pdf("PFI.pdf",height = 8,width = 16,onefile = F)
ggsurvplot(fit, data = dat, risk.table = TRUE, conf.int = F, pval=TRUE, break.time.by = 365)
dev.off()
#
dat$event=as.numeric(dat$DSS)#dat$event=as.numeric(dat$DFS)#dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$DSS.time)#dat$time=as.numeric(dat$DFS.time)#dat$time=as.numeric(dat$time)
fit <- survfit(Surv(time, event) ~cnv4, data = dat,type  = "kaplan-meier",conf.type = "log-log")
pdf("DSS.pdf",height = 8,width = 16,onefile = F)
ggsurvplot(fit, data = dat, risk.table = TRUE, conf.int = F, pval=TRUE, break.time.by = 365)
dev.off()
#
library(survival)
library(survminer)

model <- coxph(s, data =dat)#datd
summary(model,data=dat)#datd
options(scipen=1)


# 对于生存数据预后模型评价很多采用C-index ，但c-index展示没有roc曲线图来的直观
new_dat=dat
new_dat$event=as.numeric(new_dat$event)
new_dat$time=as.numeric(new_dat$time)

plot(fp)
fp <- predict(model,new_dat) ;boxplot(fp)#默认lp
as.data.frame(fp)[1,]
#分开的,type="term"
basehaz(model) 
library(Hmisc)
options(scipen=200)
with(new_dat,rcorr.cens(fp,Surv(time, event)))
# http://dni-institute.in/blogs/cox-regression-interpret-result-and-predict/
# ?surv_categorize
library(cowplot)
library(pheatmap)
  library(survival)
  library(survminer)














zh=combn(names(table(codelete$cnv4)),2)
zh
#testpn=sapply(1:ncol(zh), function(i) roc.test(roc.list[[zh[1,i]]],roc.list[[zh[2,i]]]))#全比较
testpn=lapply(1:length(zh[1,]), function(i) c(zh[1,i],zh[2,i]))#全比较
testpn=lapply(1:(length(names(table(codelete$cnv4)))-1), function(i) c(zh[1,i],zh[2,i]))#全比较1
i=1
d

list(c("No changes","Shallow Deletion"),c("No changes","MIR31HG: Deep Deletion"),c("No changes","MIR31HG: Shallow Deletion"),c("No changes","MIR31HG: Diploid"),c("No changes","MIR31HG: Gain"),c("No changes","MIR31HG: Amplification"))