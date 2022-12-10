# (RNA-seq OR "scRNA"OR "snRNA") AND (((("immune checkpoint") OR (pd-1)) OR (pd-L1)) OR (CTLA-4)) "PD" "SD" "PR" "CR" intitle:pancreatic
#还没改呢 20210309这个地方要改的很多啊。。。
#quantile normalize
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
# gset=AnnoProbe::geoChina('GSE179351')
# gset=GEOquery::getGEO('GSE179351')
load("./GSE179351_eSet.Rdata")#log不平，但像是弄过了
# check the ExpressionSet
eSet=gset[[1]]#eSet=gset[[2]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)#原矩阵,途经过滤1次
head(probes_expr[,1:4])
# probes_expr=log2(probes_expr+1)
# boxplot(probes_expr,las=2)
## pheno info
if(T){
  library(data.table)
  a1=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\icipanc\GSE179351\surv.xlsx)",sheet=1)
  a2=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\icipanc\GSE179351\area.xlsx)",sheet=1)
  # a3=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\icipanc\GSE179351\TMB.xlsx)",sheet=1)
  aa=merge(a1,a2,by="ID",all=T)
  #b1b2 不纳入因为不是重点
  b1=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\icipanc\GSE179351\suppletable.xlsx)",sheet=1,skip = 1)
  b1=b1[,-c(2,3,4,5)]
  # 
  # b2=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\icipanc\GSE179351\suppletable.xlsx)",sheet=2,skip = 1)#
  
  b3=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\icipanc\GSE179351\suppletable.xlsx)",sheet=3,skip = 1)#
  b3=unite(b3,"title",c("Patient ID","Biopsy"),sep = "-")
  b3$title=toupper(b3$title)
  b3$title=gsub(" ","",b3$title)
  b4=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\icipanc\GSE179351\suppletable.xlsx)",sheet=4)
  colnames(b4)=toupper(colnames(b4))
  b5=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\icipanc\GSE179351\suppletable.xlsx)",sheet=5)
  colnames(b5)=toupper(colnames(b5))
  bb=plyr::rbind.fill(b4,b5)
}
phenoDat <- pData(eSet)
phenoDat$title=toupper(phenoDat$title)
phenoDat$title=gsub(" ","",phenoDat$title)
phenoDat$ID=stringr::str_split(phenoDat$title,"-",simplify = T)[,1]
phenoDat$ID=as.numeric(phenoDat$ID)
phenoDat=merge(aa,phenoDat,by="ID",all=T)
phenoDat=merge(phenoDat,b3,by="title",all=T)
phenoDat[is.na(phenoDat$ID),]$ID=stringr::str_split(phenoDat[is.na(phenoDat$ID),]$title,"-",simplify = T)[,1]
phenoDat=merge(phenoDat,bb,by="ID",all=T)
phenoDat=merge(phenoDat,b1,by="ID",all=T)
phenoDat$`Out of Field Lesion Best Response RECIST`=ifelse(is.na(phenoDat$`Out of Field Lesion Best Response RECIST`),phenoDat$Response,phenoDat$`Out of Field Lesion Best Response RECIST`)
phenoDat$`Out of Field Lesion Best Response RECIST`=ifelse(is.na(phenoDat$`Out of Field Lesion Best Response RECIST`),phenoDat$`Reason for Coming Off Treatment Prior to Radiation(NE because of toxic)`,phenoDat$`Out of Field Lesion Best Response RECIST`)
phenoDat$`Out of Field Lesion Best Response RECIST`=ifelse(is.na(phenoDat$`Out of Field Lesion Best Response RECIST`),"NE",phenoDat$`Out of Field Lesion Best Response RECIST`)
phenoDat$Response=phenoDat$`Out of Field Lesion Best Response RECIST`
phenoDat$`Out of Field Lesion Best Response RECIST`=NULL
phenoDat$`Reason for Coming Off Treatment Prior to Radiation(NE because of toxic)`=NULL

# table(is.na(phenoDat$Response))
# table(is.na(phenoDat$`Out of Field Lesion Best Response RECIST`))

phenoDat$OS.time=phenoDat$`OS (Months)`*30
phenoDat$OS=phenoDat$`OS Censor`
phenoDat$PFS=phenoDat$`PFS Censor`
phenoDat$PFS.time=phenoDat$`PFS (Months)`*30
phenoDat$sample=phenoDat$title
duplicated(na.omit(phenoDat$geo_accession))
openxlsx::write.xlsx(phenoDat,file = "phenoDat_GSE179351.xlsx")

expr_raw=fread(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\icipanc\GSE179351\GSE179351_RawCountsForAllSamples.txt.gz)")

expr_tpm=expr_raw
expr_tpm$mean=rowMeans(expr_tpm[,-c(1:11)])
# expr_tpm=subset(expr_tpm,expr_tpm$mean!=0)
table(expr_tpm$gene.type)
expr_tpm=subset(expr_tpm,expr_tpm$chr%in%c(1:22,"X","Y"))
table(expr_tpm$chr)
expr_tpm=dplyr::arrange(expr_tpm,desc(mean))
expr_tpm=expr_tpm[!duplicated(expr_tpm$ensembl.gene.id),]
expr_tpm=expr_tpm[!duplicated(expr_tpm$gene.symbol),]
expr_tpm$mean=NULL
expr_tpm_realtion=expr_tpm[,1:11]
expr_tpm=expr_tpm[,-c(1:11)]
colnames(expr_tpm)=gsub("\\.","-",colnames(expr_tpm))
colnames(expr_tpm)=toupper(colnames(expr_tpm))
colnames(expr_tpm)=str_sub(colnames(expr_tpm),2, -1)
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
expr_tpm <- apply(expr_tpm,2,function(x) countToTpm(x,effLen=expr_tpm_realtion$length))
expr_tpm_ensg=expr_tpm
expr_tpm_ensg=as.data.frame(expr_tpm_ensg)
rownames(expr_tpm_ensg)=expr_tpm_realtion$ensembl.gene.id
expr_tpm_symbol=expr_tpm
expr_tpm_symbol=as.data.frame(expr_tpm_symbol)
rownames(expr_tpm_symbol)=expr_tpm_realtion$gene.symbol

expr_raw$mean=rowMeans(expr_raw[,-c(1:11)])
expr_raw=subset(expr_raw,expr_raw$mean!=0)
table(expr_raw$gene.type)
expr_raw=subset(expr_raw,expr_raw$chr%in%c(1:22,"X","Y"))
table(expr_raw$chr)
expr_raw=dplyr::arrange(expr_raw,desc(mean))
expr_raw=expr_raw[!duplicated(expr_raw$ensembl.gene.id),]
expr_raw=expr_raw[!duplicated(expr_raw$gene.symbol),]
expr_raw$mean=NULL
expr_raw_realtion=expr_raw[,1:11]
expr_raw=expr_raw[,-c(1:11)]
colnames(expr_raw)=gsub("\\.","-",colnames(expr_raw))
colnames(expr_raw)=toupper(colnames(expr_raw))
colnames(expr_raw)=str_sub(colnames(expr_raw),2, -1)
expr_raw_ensg=expr_raw
expr_raw_ensg=as.data.frame(expr_raw_ensg)
rownames(expr_raw_ensg)=expr_raw_realtion$ensembl.gene.id
expr_raw_symbol=expr_raw
expr_raw_symbol=as.data.frame(expr_raw_symbol)
rownames(expr_raw_symbol)=expr_raw_realtion$gene.symbol

expr_DESeq2_NormalizedCounts=fread(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\icipanc\GSE179351\GSE179351_DESeq2_NormalizedCountsForAllSamples.txt.gz)")

expr_DESeq2_NormalizedCounts$mean=rowMeans(expr_DESeq2_NormalizedCounts[,-c(1:11)])
expr_DESeq2_NormalizedCounts=subset(expr_DESeq2_NormalizedCounts,expr_DESeq2_NormalizedCounts$mean!=0)
table(expr_DESeq2_NormalizedCounts$gene.type)
expr_DESeq2_NormalizedCounts=subset(expr_DESeq2_NormalizedCounts,expr_DESeq2_NormalizedCounts$chr%in%c(1:22,"X","Y"))
table(expr_DESeq2_NormalizedCounts$chr)
expr_DESeq2_NormalizedCounts=dplyr::arrange(expr_DESeq2_NormalizedCounts,desc(mean))
expr_DESeq2_NormalizedCounts=expr_DESeq2_NormalizedCounts[!duplicated(expr_DESeq2_NormalizedCounts$ensembl.gene.id),]
expr_DESeq2_NormalizedCounts=expr_DESeq2_NormalizedCounts[!duplicated(expr_DESeq2_NormalizedCounts$gene.symbol),]
expr_DESeq2_NormalizedCounts$mean=NULL
expr_DESeq2_NormalizedCounts_realtion=expr_DESeq2_NormalizedCounts[,1:11]
expr_DESeq2_NormalizedCounts=expr_DESeq2_NormalizedCounts[,-c(1:11)]
colnames(expr_DESeq2_NormalizedCounts)=gsub("\\.","-",colnames(expr_DESeq2_NormalizedCounts))
colnames(expr_DESeq2_NormalizedCounts)=toupper(colnames(expr_DESeq2_NormalizedCounts))
colnames(expr_DESeq2_NormalizedCounts)=str_sub(colnames(expr_DESeq2_NormalizedCounts),2, -1)
expr_DESeq2_NormalizedCounts_ensg=expr_DESeq2_NormalizedCounts
expr_DESeq2_NormalizedCounts_ensg=as.data.frame(expr_DESeq2_NormalizedCounts_ensg)
rownames(expr_DESeq2_NormalizedCounts_ensg)=expr_DESeq2_NormalizedCounts_realtion$ensembl.gene.id
expr_DESeq2_NormalizedCounts_symbol=expr_DESeq2_NormalizedCounts
expr_DESeq2_NormalizedCounts_symbol=as.data.frame(expr_DESeq2_NormalizedCounts_symbol)
rownames(expr_DESeq2_NormalizedCounts_symbol)=expr_DESeq2_NormalizedCounts_realtion$gene.symbol


c(library(survival),library(survminer),library(ggtext),library(cowplot))
library(Hmisc)

library(cowplot)
library(pheatmap)
# after_multiple_surviaval(genes="ENSG00000140718",rna.expr=TCGA_tpm_ensg,metaMatrix=pcawg_clinical_CDR,i) # returns list
genes="ENSG00000140718"
# rna.expr=expr_tpm_ensg
rna.expr=expr_raw_ensg
# rna.expr=expr_DESeq2_NormalizedCounts_ensg
rna.expr=log2(rna.expr+1)
rna.expr=limma::normalizeBetweenArrays(rna.expr)
rna.expr=as.data.frame(rna.expr)
phenoDat_exp=subset(phenoDat,!is.na(phenoDat$geo_accession))
phenoDat_exp$therapy=as.factor(stringr::str_split(phenoDat_exp$sample,"-",simplify = T,n = 2)[,2])
# phenoDat_exp=subset(phenoDat_exp,therapy%in%c("PRE-TX","PRE-XRT"))
# phenoDat_exp=dplyr::arrange(phenoDat_expw,phenoDat_expw$therapy)
# phenoDat_exp=phenoDat_exp[!duplicated(phenoDat_exp$ID),]
# phenoDat_exp=subset(phenoDat_exp,!is.na(phenoDat_exp$KRAS))
# phenoDat_exp=subset(phenoDat_exp,!is.na(phenoDat_exp$`Out of Field Lesion Best Response RECIST`))
# phenoDat_exp=subset(phenoDat_exp,phenoDat_exp$`cancer type:ch1`=="PDAC")
phenoDat_exp_tx=subset(phenoDat_exp,grepl("tx",phenoDat_exp$title,ignore.case = T))
# phenoDat_exp_tx=subset(phenoDat_exp_tx,OS.time>60)
# phenoDat_exp_tx=subset(phenoDat_exp_tx,PFS.time>30)
# phenoDat_exp_prexrt=subset(phenoDat_exp,grepl("pre-xrt",phenoDat_exp$title,ignore.case = T))
# phenoDat_exp_postxrt=subset(phenoDat_exp,grepl("post",phenoDat_exp$title,ignore.case = T))

metaMatrix=phenoDat_exp_tx
dat=before_multiple_surviaval(genes, rna.expr, metaMatrix)
# metaMatrix=phenoDat_exp_postxrt
# metaMatrix=phenoDat_exp_prexrt


dat$FTO=ifelse(dat$ENSG00000140718>median(dat$ENSG00000140718,na.rm = T),"FTO-high","FTO-low")
d=dat
d=subset(d,d$Response!="NE")
d$Response=factor(d$Response,levels = c("CR","PR","SD","PD"),labels = c("R","R","R","PD"),ordered = F)#,"NE"
#####

library(tidyverse)
library(RColorBrewer)
set.seed(123)

theme_big_simple1 <- function(){ 
  theme_bw(base_size=16, base_family="") %+replace% 
    theme(
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      legend.title=element_text(vjust = 1,face="bold",family="sans",size=24),
      legend.text=element_text(vjust = 1,face="bold",family="sans",size=20),
      legend.justification = "right",
      legend.key.width = unit(40, "pt"),
      axis.line = element_line(color = "black", size = 1, linetype = "solid"),
      axis.ticks = element_line(colour = "black", size = 1),
      panel.background  = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),#panel.grid.major = element_line(color="lightgrey"),
      panel.border = element_blank(),
      legend.position = ("bottom"),
      plot.title = element_text(face="bold",family="sans",size=24, hjust = 0.0, vjust = 1.75),
      axis.text.x = element_text(face="bold",family="sans",color="black", size=20, margin = margin(t = 4, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(face="bold",family="sans",hjust = 1,color="black", size=20, margin = margin(t = 0, r = 4, b = 0, l = 0)),
      axis.title.x = element_text(face="bold",family="sans",margin = margin(t = 10, r = 0, b = 0, l = 0), angle = 0, size = 24),
      axis.title.y = element_text(face="bold",family="sans",hjust = 1,margin = margin(t = 0, r = 10, b = 0, l = 0), angle = 90, size = 24),
      axis.ticks.length = unit(0.20, "cm"),
      strip.background = element_rect(color="black", size=1, linetype="solid"),
      strip.text.x = element_text(face="bold",family="sans",size = 20, color = "black"),
      strip.text.y = element_text(face="bold",family="sans",size = 20, color = "black")
    )}



set.seed(123)
library(ggstatsplot)

# association test (or contingency table analysis)
g=ggbarstats(d, x = Response, y =FTO, results.subtitle = F,label="both",label.args = list(alpha = 0,label.size=NA))+theme_big_simple1()+
  labs(x="", y="")+
  ggsci::scale_fill_nejm()
if(F){ 
grouped_ggbarstats(d, x = Response, y =FTO,grouping.var = `cancer type:ch1`)

g <- ggplot(d, aes(FTO)) + geom_bar(aes(fill = Response),position='fill')+#position = position_stack(reverse = TRUE)) +
  theme(legend.position = "top")+  labs(x="", y=" ")+
  ggsci::scale_fill_nejm()+
  theme_classic()+
  theme(legend.position = "top",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1))+theme_big_simple1()
}
pdf(file= paste0("Figure.response.pdf"),width = 4.2,height=8)
g
dev.off()
ggsave(filename = paste0("Figure.response.pdf"),
       print(g),
       device = 'pdf',
       dpi = 1200)
# pvalue <- chisq.test(c(125,24,105,44,ncol=2))$p.value #卡方检验
# library(plyr)
# ggplot(a,aes(Var1,percent,fill=Var2))+
#   geom_bar(stat="identity",position = position_stack())+
#   scale_fill_manual(values = c("#DB423E","#008ECA"),label=c("SD/PD","CR/PR"))+
#   scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
#   labs(x="Risk",y="Percent Weidght",
#        fill="")+
#   geom_text(aes(label=label),vjust=3,size=6,color="black")+
#   annotate(geom = "text",
#            cex=6,
#            x=1.5, y=105, # 根据自己的数据调节p value的位置
#            label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
#            color="black")+
#   theme_classic()+
#   theme(legend.position = "top",
#         legend.text = element_text(size=12),
#         axis.text = element_text(size=12),
#         axis.title = element_text(size=12))



# multiple_surviaval(genes,dat,name = "FTO")
# multiple_surviaval(x,dat)
x="ENSG00000140718"

before_multiple_surviaval <- function(genes, rna.expr, metaMatrix){
  library(survminer)
  library(survival)
  samples = intersect(colnames(rna.expr), metaMatrix$sample)
  exprDa=rna.expr[genes,samples]
  clinicalDa=metaMatrix[match(samples,metaMatrix$sample),]
  colnames(exprDa)==rownames(clinicalDa)
  dat=cbind(t(exprDa),clinicalDa)
  return(dat)
}

multiple_surviaval <- function(x,dat,name){
  dir.create("./piliangshengcuntu/")
  dir.create(paste0("./piliangshengcuntu/",name))
  dir.create(paste0("./piliangshengcuntu/",name,"/mediancut"))
  dir.create(paste0("./piliangshengcuntu/",name,"/bestcut"))
  
  dat$event=as.numeric(dat$OS)
  dat$time=as.numeric(dat$OS.time)
  # dat$event=as.numeric(dat$PFS)
  # dat$time=as.numeric(dat$PFS.time)
  dat=dat[!is.na(dat$event),]
  dat=dat[!is.na(dat$time),]
  # dat=dat[dat$time>60,]
  # dat=dat[rownames(dat)!="12-PRE-TX",]#也不算outlier
  
  if(F){
  s=as.formula(paste0("survival::Surv(time, event) ~",as.character(x)))
  model <- survival::coxph(formula=s, data =dat,nocenter = T)
  summary(model,data=dat)
  # options(scipen=1)
  new_dat=dat
  # fp <- predict(model,new_dat)
  risk=new_dat
  risk$fp=fp
}
  library(cowplot)
  library(pheatmap)
  new_dat=dat
  risk=new_dat
  risk$fp=dat[[x]]
  res.cut.va <- surv_cutpoint(risk, time = "time", event = "event",minprop=0.3,
                              variables =c("fp"))
  summary(res.cut.va)
  res.cat.va <- survminer::surv_categorize(res.cut.va)
  fit <- survival::survfit(survival::Surv(time, event) ~fp, data = res.cat.va)
  
  # pdf(file=paste0("./piliangshengcuntu/bestcut","Figure.",x,"_%03d","_best_survplot.pdf"),height = 6,width = 6,onefile = T)
  edr1=survminer::ggsurvplot(fit,
                             data = res.cat.va,
                             # legend.title = paste0("OS","-",name,"-",ensg_symbol_trans[ensg_symbol_trans$id==x,"gene_name"]),#as.character(choose_gene[choose_gene$Gene%in% as.character(x),"GeneName"]),
                             legend.labs = c( "High","Low"),
                             risk.table = T,#循环时不能用否则会报错，单个人可以T
                             conf.int = F,
                             xlim=c(0,540),
                             break.x.by=90,
                             pval=TRUE,
                             tables.height = 0.2,
                             tables.theme = theme_cleantable(),
                             palette = c("#E7B800","#2E9FDF"),
                             ggtheme = theme_bw()
  )
  g=edr1
  g
  # dev.off()
  if(F){
    ggsurv=edr1
    g2 <- ggplotGrob(ggsurv$plot)
    g3 <- ggplotGrob(ggsurv$table)
    min_ncol <- min(ncol(g2), ncol(g3))
    g <- gridExtra::gtable_rbind(g2[, 1:min_ncol], g3[, 1:min_ncol], size="last")
    g$widths <- grid::unit.pmax(g2$widths, g3$widths)
    grid::grid.newpage()
    grid::grid.draw(g)
  }
  # ggsave(filename = paste0("Figure.",x,"_best_survplot.pdf"),
  #        plot=print(edr1, newpage = FALSE),
  #        path = paste0("./piliangshengcuntu/",stringr::str_split(multiplesurvival,"\\/",simplify=T)[,8],"/bestcut"),
  #        device = 'pdf',
  #        width = 8,height=8,dpi = 600)
  ggsave(filename = paste0("Figure.",x,"_",name,"_best_survplot.pdf"),
         print(g$plot),
         path = paste0("./piliangshengcuntu/",name,"/bestcut"),
         device = 'pdf',
         width = 4,height=4,dpi = 600)
  
  # dev.off()
  # pdf(file=paste0("./piliangshengcuntu/","Figure.",x,"_km_median_survplot.pdf"),height = 8,width = 8,onefile = F)
  res.cut.va$cutpoint$cutpoint=median(res.cut.va$data$fp)
  res.cat.va <- surv_categorize(res.cut.va)
  fit <- survfit(survival::Surv(time, event) ~fp, data = res.cat.va)
  edr2=survminer::ggsurvplot(fit,
                             data = res.cat.va,
                             legend.title = paste0("OS","-",name,"-",ensg_symbol_trans[ensg_symbol_trans$id==x,"gene_name"]),#as.character(choose_gene[choose_gene$Gene%in% as.character(x),"GeneName"]),
                             legend.labs = c( "High","Low"),
                             risk.table = F,#循环时不能用否则会报错，单个人可以T
                             conf.int = TRUE,
                             xlim=c(0,1830),
                             break.x.by=365,
                             pval=TRUE,
                             tables.height = 0.2,
                             tables.theme = theme_cleantable(),
                             palette = c("#E7B800", "#2E9FDF"),
                             ggtheme = theme_bw()
  )
  g=edr2
  
  if(F){
    ggsurv=edr2
    g2 <- ggplotGrob(ggsurv$plot )
    g3 <- ggplotGrob(ggsurv$table)
    min_ncol <- min(ncol(g2), ncol(g3))
    g <- gridExtra::gtable_rbind(g2[, 1:min_ncol], g3[, 1:min_ncol], size="last")
    g$widths <- grid::unit.pmax(g2$widths, g3$widths)
    grid::grid.newpage()
    grid::grid.draw(g)
  }
  # ggsave(filename = paste0("Figure.",x,"_km_median_survplot.pdf"),
  #        plot=print(edr2, newpage = FALSE),
  #        path = paste0("./piliangshengcuntu/",stringr::str_split(multiplesurvival,"\\/",simplify=T)[,8],"/mediancut"),
  #        device = 'pdf',
  #        width = 8,height=8,dpi = 600)
  ggsave(filename = paste0("Figure.",x,"_",name,"_km_median_survplot.pdf"),
         print(g$plot),
         path = paste0("./piliangshengcuntu/",name,"/mediancut"),
         device = 'pdf',
         width = 4,height=4,dpi = 600)
  # dev.off()
  # sprintf("p%d.png", i)
}

after_multiple_surviaval <- function(genes, rna.expr, metaMatrix,i){
  dat=before_multiple_surviaval(genes, rna.expr[[i]], metaMatrix)
  multiple_surviaval(genes,dat,names(rna.expr)[[i]])
}

probesetvar = apply(expr_raw_ensg, 1, var)  #表达变化大的基因
ord = order(probesetvar, decreasing=TRUE)[1:1000]   #前200个基因，或者更多
pca = prcomp(t(expr_tpm_ensg[ord,]), scale=TRUE)
ss=summary(pca)
library(factoextra)
library(FactoMineR)
pdf(file='pca ici pdac.pdf',height = 10,width = 10)#var,quali,all,none
fviz_pca_ind(pca, #label="none", 
             habillage=as.factor(stringr::str_split(colnames(expr_tpm_ensg),"-",simplify = T,n = 2)[,2]),
             addEllipses=TRUE,labelsize=4,pointsize=2, ellipse.level=0.95, palette = NULL)+ggtitle("")+theme(text = element_text(size = 25),legend.position = "top") #+ coord_fixed(1)
dev.off()


# table(expr_raw$gene.symbol==expr_DESeq2_NormalizedCounts$gene.symbol)#TRUE 76780 
# expr_raw_realtion=expr_raw[,1:11]
# expr_raw_ensg=expr_raw[,-c(1:11)]
# rownames(expr_raw_ensg)=expr_raw$ensembl.gene.id
# table(duplicated(expr_raw$ensembl.gene.id))
# table(duplicated(expr_raw$gene.symbol))
# expr_raw_symbol=expr_raw[,-c(1:11)]
# expr_DESeq2_NormalizedCounts=expr_DESeq2_NormalizedCounts[,-c(1:11)]

if(eSet@annotation=="GPL570"){
  system("tar -xvf ./GSE179351/GSE179351_RAW.tar -C ./GSE179351/")
  celFiles <- list.celfiles('./GSE179351/',full.name=TRUE,listGzipped = T)
  affyRaw <- read.celfiles(celFiles)
  # 提取矩阵并做normalization 
  #rma 自带 背景校正方法：rma，标准化方法：quantile，汇总方法：medianpolish
  #不需要重复quantile，当然limma要干嘛随便
  eset <- rma(affyRaw)#不存在直接可以复现的quantile，因为这是用C语言写的
  # 检查数据
  probes_expr_rma=exprs(eset)
  dim(probes_expr_rma)
  #
  pdf("GSE179351_boxplot.pdf",width=100)
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
# phenoDat_outcome=read.csv("./GSE179351/GSE179351_outcome.csv",header = T)
# phenoDat_outcome=as.data.frame(phenoDat_outcome)
#
# phenoDat=merge(phenoDat,phenoDat_outcome,by.x="patient:ch1",by.y="Tumor.ID")
# rownames(phenoDat)=phenoDat$geo_accession
## check GPL and annotate the probes to genes.
if(F){
  (gpl=eSet@annotation)
  checkGPL(gpl)
  printGPLInfo(gpl)
  probe2gene=idmap(gpl=gpl,type = "soft")#soft的合并性差落后，建议首选bioc,但如果结合了更新就不一样了,例如ta
  head(probe2gene)
}
#为了图方便我用了直接下载的矩阵!!!!!
probe2gene=as.data.frame(probes_expr)
probe2gene$symbol=rownames(probes_expr)
probe2gene$ID=rownames(probes_expr)
#可能不行啊（symbol滤过了少了很多啊）
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
# 修剪+median !!!!!!!!!!!!!!!!!!!!!!!!
dim(phenoDat)
# phenoDat$time=as.numeric(phenoDat$time)
# phenoDat=phenoDat[phenoDat$time>=1,]
# dim(phenoDat)
phenoDat=phenoDat[phenoDat$characteristics_ch2.1=="tissue type: Primary",]
phenoDat=phenoDat[!is.na(phenoDat$`death_event_1death_0censor:ch2`),]

phenoDat$OS=as.numeric(phenoDat$`death_event_1death_0censor:ch2`)
phenoDat$OS.time=as.numeric(phenoDat$`survival_months:ch2`)
phenoDat_GSE179351_rna_125_keepall=phenoDat

genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
genes_expr_mean_GSE179351_rna_125=t(scale(t(genes_expr_with_clinical)))

phenoDat=subset(phenoDat,OS.time>=1)

genes_expr_with_clinical=genes_expr[,which(colnames(genes_expr) %in% phenoDat$geo_accession)]
dim(genes_expr_with_clinical)
genes_expr_GSE179351_rna=genes_expr_with_clinical
phenoDat_GSE179351_rna=phenoDat
probes_expr_GSE179351_rna=probes_expr
genes_expr_mad_GSE179351_rna=t(RobScale(t(genes_expr_with_clinical)))
genes_expr_mean_GSE179351_rna=t(scale(t(genes_expr_with_clinical)))


#我觉得不用删掉，只是缺东西而已。全删除就是97个、


#
save(phenoDat_GSE179351_rna_125_keepall,genes_expr_mean_GSE179351_rna_125,phenoDat_GSE179351_rna,probes_expr_GSE179351_rna,genes_expr_mad_GSE179351_rna,genes_expr_mean_GSE179351_rna,genes_expr_GSE179351_rna,
     file="GSE179351_after_soft.Rdata")
#(file="GSE179351_after_bioc.Rdata"))
write.csv(phenoDat_GSE179351_rna_125_keepall,file="phenoDat_GSE179351_rna_125_keepall.csv",quote = T)
write.csv(phenoDat_GSE179351_rna,file="phenoDat_GSE179351_rna.csv",quote = T)
pdf("GSE179351_mad.pdf",width = 100)
par(mfrow = c(2,1));
boxplot(probes_expr)
genes_expr_mean=t(scale(t(genes_expr)))
boxplot(genes_expr_mean)
dev.off()

