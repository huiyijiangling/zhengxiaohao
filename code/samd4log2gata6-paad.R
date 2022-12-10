#20210713 v2 mir31hg改为二倍体型
rm(list=ls())
options(stringsAsFactors = F)
gc()
# library(GDCRNATools)

SMAD4txt=read.table(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/hic scripts/samd_gata6/paad-SMAD4-GATA6.txt",sep = "\t",header = T)
colnames(SMAD4txt)
codelete=SMAD4txt
stringr::str_split(codelete$Copy.Number.Alterations,": ",simplify = T)
codelete$GATA6status=stringr::str_split(codelete$Copy.Number.Alterations,": ",simplify = T)[,2]
# MIR31HGtxt=read.csv(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/hic scripts/tcga 4 gene/mir31hg.txt",sep = "\t",header = T)
# codelete=merge(codelete,MIR31HGtxt,by="Sample.Id")
# codelete$cnv=ifelse(codelete$cnv=="Diploid Diploid Diploid","No changes",
#                     ifelse(codelete$cnv=="Shallow Deletion Shallow Deletion Shallow Deletion","Shallow Deletion","Deep Deletion"))
# table(codelete$Copy.Number.Alterations)
# table(codelete$MIR31HG..Putative.copy.number.alterations.from.GISTIC)
# codelete=subset(codelete,codelete$Copy.Number.Alterations=="MIR31HG: Diploid")
# codelete=subset(codelete,!codelete$cnv=="Shallow Deletion")#不要shallow del
# codelete=subset(codelete,!(cnv %in% "No changes"&MIR31HG..Putative.copy.number.alterations.from.GISTIC != "Diploid"))#9
codelete$cnv=codelete$SMAD4..Putative.copy.number.alterations.from.GISTIC#ifelse(codelete$cnv!="Deep Deletion",codelete$cnv,codelete$Copy.Number.Alterations)
codelete=subset(codelete,cnv%in%c("Deep Deletion","Diploid"))
codelete$cnv3=codelete$cnv#ifelse(codelete$cnv!="Deep Deletion",codelete$cnv,codelete$Copy.Number.Alterations)
box=codelete
if(max(box$GATA6..mRNA.Expression..RSEM..Batch.normalized.from.Illumina.HiSeq_RNASeqV2.)>50){
  box$number1=as.numeric(log2(box$GATA6..mRNA.Expression..RSEM..Batch.normalized.from.Illumina.HiSeq_RNASeqV2.+1))
}
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
write.csv(box,file="PAAD-SMAD4-GATA6.csv",quote=T)

shapiro.test#3-5000
lillie.test#>5000
library(nortest)
table(box_surv$cnv3)
tapply(box$number1,box$cnv3,lillie.test)#不正态
bartlett.test(ENSG00000171889~cnv,data=box)#方差齐性

pdf('Figure A.Boxplot of MIR31HG baed on status of codelete samd4log2gata6-paad.pdf',width = 8,height=8)
if(require('ggpubr')){
  library(ggpubr)
  # google search : ggpubr boxplot add p-value
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  p <- ggboxplot(box, x = "GATA6status", y = "number1",fill="cnv3",
                 font.label = list(size = 1, color = "black"),outlier.shape = NA#,palette=c(1,2,3),order=c("No changes","Shallow Deletion","Deep Deletion")#"#00AFBB","#E7B800",
  )+#add = "jitter"
    # scale_color_manual(values=c(1,2,3,4))  +#"#E69F00", "#56B4E9"
    coord_cartesian(ylim = c(0.0000, 20))+
    # theme(legend.position="none")+
    guides(fill = guide_legend(title = 'SMAD4'))+
    labs(title="",x="GATA6", y = "Expression")
  #+#tag="ENSG00000000003"+
  
  # ?geom_dotplot
  # geom_dotplot(binaxis='y', stackdir='center', stackratio=0.5, dotsize=0.2,fill = c("#999999"))
  #  Add p-value
  p + stat_compare_means(,method = "kruskal.test", label.y = 20 )#aes(group=cnv3)
  # stat_compare_means(aes(group=cnv3),method = "kruskal.test", label.y = as.numeric(seq(14.7,200,by=0.7)))#+ #备注3
  
  # ,comparisons =combn(1:5, 2, FUN = list)
  #   stat_compare_means(comparisons = list(c("No changes","Shallow Deletion"),c("No changes","Deep Deletion"),c("Shallow Deletion","Deep Deletion")),method = "wilcox.test",label.y = as.numeric(seq(100,850,by=20)))+ #备注3
  #   stat_summary(aes(x = cnv3, y = number1),fun= "mean", geom = "point",position = position_dodge(0.75), shape = 23, size = 1, fill = "pink")
  # 
}
dev.off()
testpn
write.csv(box,file="co-CDKN2A-CDKN2B-MTAP-tcga-cnv3.csv",quote=T)

cellall=readxl::read_xlsx("C:/Users/zxh/Desktop/R/original cellPAAD.xlsx",sheet=3)
table(duplicated(substr(box$Sample.Id,1,12)))#false
box$bcr_patient_barcode=substr(box$Sample.Id,1,12)
box_surv=merge(box,cellall,by.x="bcr_patient_barcode",by.y = "bcr_patient_barcode")
box_surv=subset(box_surv,!is.na(box_surv$OS.time))
box_surv=subset(box_surv,box_surv$OS.time>30)
dat=box_surv
library(survival)
library(survminer)

box_surv=subset(box_surv,!ajcc_pathologic_tumor_stage %in% c("[Discrepancy]","[Not Applicable]","[Not Available]","I/II NOS","IS","Stage 0","Stage X"))
box_surv$stage=factor(box_surv$ajcc_pathologic_tumor_stage,c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB","Stage IIC","Stage III","Stage IIIA","Stage IIIB","Stage IIIC","Stage IV","Stage IVA","Stage IVB","Stage IVC"), labels = c("Stage I","Stage I","Stage I","Stage II","Stage II","Stage II","Stage II","Stage III","Stage III","Stage III","Stage III","Stage IV","Stage IV","Stage IV","Stage IV"),ordered=T)

dat=matchdata
dat$event=as.numeric(dat$OS)#dat$event=as.numeric(dat$DFS)#dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$OS.time)#dat$time=as.numeric(dat$DFS.time)#dat$time=as.numeric(dat$time)
#
fit <- survfit(Surv(time, event) ~cnv3, data = dat,type  = "kaplan-meier",conf.type = "log-log")
pdf("os.pdf",height = 8,width = 16,onefile = F)
ggsurvplot(fit, data = dat, risk.table = TRUE, conf.int = F, pval=TRUE, break.time.by = 365)
dev.off()
#
dat$event=as.numeric(dat$DFI)#dat$event=as.numeric(dat$DFS)#dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$DFI.time)#dat$time=as.numeric(dat$DFS.time)#dat$time=as.numeric(dat$time)
fit <- survfit(Surv(time, event) ~cnv3, data = dat,type  = "kaplan-meier",conf.type = "log-log")
pdf("dfi.pdf",height = 8,width = 16,onefile = F)
ggsurvplot(fit, data = dat, risk.table = TRUE, conf.int = F, pval=TRUE, break.time.by = 365)
dev.off()
#
dat$event=as.numeric(dat$PFI)#dat$event=as.numeric(dat$DFS)#dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$PFI.time)#dat$time=as.numeric(dat$DFS.time)#dat$time=as.numeric(dat$time)
fit <- survfit(Surv(time, event) ~cnv3, data = dat,type  = "kaplan-meier",conf.type = "log-log")
pdf("PFI.pdf",height = 8,width = 16,onefile = F)
ggsurvplot(fit, data = dat, risk.table = TRUE, conf.int = F, pval=TRUE, break.time.by = 365)
dev.off()
#
dat$event=as.numeric(dat$DSS)#dat$event=as.numeric(dat$DFS)#dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$DSS.time)#dat$time=as.numeric(dat$DFS.time)#dat$time=as.numeric(dat$time)
fit <- survfit(Surv(time, event) ~cnv3, data = dat,type  = "kaplan-meier",conf.type = "log-log")
pdf("DSS.pdf",height = 8,width = 16,onefile = F)
ggsurvplot(fit, data = dat, risk.table = TRUE, conf.int = F, pval=TRUE, break.time.by = 365)
dev.off()
#psm 大型社死现场
box_surv$type
paste(names(table(box_surv$ajcc_pathologic_tumor_stage)),collapse = '","')

preBL <- CreateTableOne(vars=c("cnv3","age_at_initial_pathologic_diagnosis","type","gender","race","histological_grade","stage"),
                        strata="cnv3",data=box_surv,includeNA=T,
                        factorVars=c("cnv3","type","gender","race","histological_grade","stage"))
# treat是感兴趣变量,re78为结局变量
print(preBL,showAllLevels = TRUE)

table1 <- print(preBL, 
                printToggle = FALSE, 
                noSpaces = TRUE)
kable(table1,  
      align = 'c', 
      caption = 'Table 1: Comparison of unmatched samples')
box_surv222222=subset(box_surv,is.na(box_surv$age_at_initial_pathologic_diagnosis))
age_at_initial_pathologic_diagnosis
table(box_surv$age_at_initial_pathologic_diagnosis,useNA = "ifany")
box_surv$treat=ifelse(box_surv$cnv=="Deep Deletion",1,0)
f=matchit(treat~stage+type,data=box_surv,method="nearest",exact=~type,ratio = 1,caliper=0.05)
# f=matchit(treat~age_at_initial_pathologic_diagnosis+gender+race+stage+type,data=box_surv,method="nearest",exact=~stage+type,ratio = 1,caliper=0.1)
# treat是感兴趣变量,re78为结局变量

summary(f)
# ...
# Sample Sizes:
#   Control Treated
# All           429     185
# Matched       185     185
# Unmatched     244       0
# Discarded       0       0

matchdata=match.data(f)

mBL <- CreateTableOne(vars=c("cnv3","age_at_initial_pathologic_diagnosis","type","gender","race","histological_grade","stage"),
                      strata="cnv3",data=matchdata,includeNA=T,
                      factorVars=c("cnv3","type","gender","race","histological_grade","stage"))
print(mBL,showAllLevels = TRUE)
plot(f, type = 'jitter', interactive = FALSE)
set.seed(1234)

library(foreign)
matchdata$id<-1:nrow(matchdata)
write.csv(matchdata,"matchdata.csv")

#

# model <- coxph(s, data =dat)#datd
# summary(model,data=dat)#datd
# options(scipen=1)
# 
# 
# # 对于生存数据预后模型评价很多采用C-index ，但c-index展示没有roc曲线图来的直观
# new_dat=dat
# new_dat$event=as.numeric(new_dat$event)
# new_dat$time=as.numeric(new_dat$time)
# 
# plot(fp)
# fp <- predict(model,new_dat) ;boxplot(fp)#默认lp
# as.data.frame(fp)[1,]
# #分开的,type="term"
# basehaz(model) 
# library(Hmisc)
# options(scipen=200)
# with(new_dat,rcorr.cens(fp,Surv(time, event)))
# # http://dni-institute.in/blogs/cox-regression-interpret-result-and-predict/
# # ?surv_categorize
# library(cowplot)
# library(pheatmap)
# library(survival)
# library(survminer)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# zh=combn(names(table(codelete$cnv3)),2)
# zh
# #testpn=sapply(1:ncol(zh), function(i) roc.test(roc.list[[zh[1,i]]],roc.list[[zh[2,i]]]))#全比较
# testpn=lapply(1:length(zh[1,]), function(i) c(zh[1,i],zh[2,i]))#全比较
# testpn=lapply(1:(length(names(table(codelete$cnv3)))-1), function(i) c(zh[1,i],zh[2,i]))#全比较1
# i=1
# d
# 
# list(c("No changes","Shallow Deletion"),c("No changes","MIR31HG: Deep Deletion"),c("No changes","MIR31HG: Shallow Deletion"),c("No changes","MIR31HG: Diploid"),c("No changes","MIR31HG: Gain"),c("No changes","MIR31HG: Amplification"))
geom_signif