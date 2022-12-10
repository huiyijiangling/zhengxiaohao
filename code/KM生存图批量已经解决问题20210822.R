source("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/updateName.R")# ncbi与biomart相比多map了200，等于是，ncbi》biomart》org.hs
probe2gene=choose_gene1
probe2gene$symbol=probe2gene$GeneName
forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol)
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
probe2gene=unique(probe2gene[,c("GeneName","ENSEMBL")])#bioc
# probe2gene=unique(probe2gene[,c("ID","ENSEMBL")])#soft
# probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
library(AnnoProbe)
colnames(probe2gene)=c("probe_id","symbol")
choose_gene=probe2gene
choose_gene$Gene=probe2gene$probe_id

dat=t(subset(rnaExpr,rownames(rnaExpr)%in%choose_gene$symbol))
dat=as.data.frame(dat)
choose_gene=choose_gene[choose_gene$symbol %in% colnames(dat),]
dat$sample=rownames(dat)
dat=merge(phe,dat,by="sample")
choose_gene[!choose_gene$Gene %in% colnames(dat),]

#PCAWG
i=1
if(F){
  datasetExpr=list()
  dat=list()
  #pcawg_clinical_CDR,TCGA_tpm_symbol,TCGA_tpm_ensg
  datasetExpr[[i]]=TCGA_tpm_ensg[[i]]
  phe=pcawg_clinical_CDR
  choose_gene=ensg_symbol_trans
  dat[[i]]=t(subset(datasetExpr[[i]],rownames(datasetExpr[[i]])%in%choose_gene$id))
  dat[[i]]=as.data.frame(dat[[i]])
  choose_gene=choose_gene[choose_gene$id %in% colnames(dat[[i]]),]
  dat[[i]]$sample=rownames(dat[[i]])
  dat[[i]]=merge(phe,dat[[i]],by="sample")
  choose_gene[!choose_gene$id %in% colnames(dat),]
}
library(doParallel) #
cl <- makeCluster(33)
registerDoParallel(cl)
foreach(x=survOutput005$ensg,.packages = c(library(survival),library(survminer),library(ggtext),library(cowplot))) %dopar% multiple_surviaval(x,dat) # returns list
# choose_gene[[1]]
stopCluster(cl)
#一模一样啊，。。。。。。没区别
#####################################################################################

rm(list=ls())
options(stringsAsFactors = F)
gc()

library(readxl)
# multiplesurvival <- list.files('C:/Users/zxh/Desktop/R/hic/multiple survival/up',pattern = ".xlsx",full.name=TRUE)
multiplesurvival <- list.files('C:/Users/zxh/Desktop/R/hic/multiple survival/down',pattern = ".xlsx",full.name=TRUE)

multiplesurvival
choose_gene <- lapply(multiplesurvival, function(fl) read_excel(path=fl,sheet = 1))[[1]][,1:2]



load("C:/Users/zxh/Desktop/R/paad-tcga-gtex/放进paad/calculate survival.Rdata")




library(survival)
library(survminer)
library(ggtext)
# options(warn=0)
# paste(rownames(cc),collapse = "+")
# choose_gene=unlist(lapply(multiplesurvival,function(x) x[[1]]))

# rnaCounts["ENSG00000271447",]
# dat=t(rnaExpr)#这种对于del特别多的就不合适
dat=t(rnaExpr_notfilter)
dat=as.data.frame(dat)
dat$sample=rownames(dat)
dat=merge(phe,dat,by="sample")
choose_gene[!choose_gene$Gene %in% colnames(dat),]
#CCDC177 表达过低不建议使用
#MMP28 目前版本无合适对应
choose_gene=choose_gene[choose_gene$Gene %in% colnames(dat),]

if(F){
  #另一项 杜老师临时筛选chr18上基因任务
  library(readxl)
  multiplesurvival <- list.files('C:/Users/zxh/Desktop/R/paad-tcga-gtex/放进paad/chr18',pattern = ".xlsx",full.name=TRUE)
  
  multiplesurvival
  choose_gene1 <- lapply(multiplesurvival, function(fl) read_excel(path=fl,sheet = 1))[[1]][,4]
  #单基因
  if(F){
    choose_gene1="MIR31HG"
    choose_gene1=as.data.frame(choose_gene1)
    choose_gene1$GeneName=choose_gene1$choose_gene1
  }
  choose_gene=GDCRNATools::gdcDEAnalysis(counts     = rnaCounts, 
                                         group      = phe$vital_status, 
                                         comparison = 'Dead-Alive', 
                                         method     = 'limma',#edgeR
                                         filter =F)#默认是true，不筛选则选择false
  choose_gene=choose_gene[choose_gene$symbol %in% choose_gene1$GeneName,]
  table(choose_gene$group)
  choose_gene=subset(choose_gene,group!="ncRNA" )
  choose_gene=choose_gene[rownames(choose_gene) %in% colnames(dat),]
  choose_gene$Gene=rownames(choose_gene)
  colnames()[1]="GeneName"
}

x="ENSG00000038427"
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
  dat=dat[!is.na(dat$event),]
  dat=dat[!is.na(dat$time),]
  dat=dat[which(dat$time>30),]
  # dat=dat[which(dat$PFI.time>30),]
  # if(F){
  #   s=as.formula(paste0("survival::Surv(time, event) ~",as.character(x)))
  #   model <- survival::coxph(formula=s, data =dat)
  #   # summary(model,data=dat)
  #   # options(scipen=1)
  #   new_dat=dat
  #   fp <- predict(model,new_dat)
  #   library(cowplot)
  #   library(pheatmap)
  #   risk=new_dat
  #   risk$fp=fp
  # }
  risk=dat
  risk$fp=dat[[x]]
  res.cut.va <- surv_cutpoint(risk, time = "time", event = "event",minprop=0.2,
                              variables =c("fp"))
  summary(res.cut.va)
  res.cat.va <- survminer::surv_categorize(res.cut.va)
  fit <- survival::survfit(survival::Surv(time, event) ~fp, data = res.cat.va)
  
  # pdf(file=paste0("./piliangshengcuntu/bestcut","Figure.",x,"_%03d","_best_survplot.pdf"),height = 6,width = 6,onefile = T)
  edr1=survminer::ggsurvplot(fit,
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
                             palette = c("#E7B800","#2E9FDF"),
                             ggtheme = theme_bw()
  )
  g=edr1
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
# ENSG00000140718 FTO
# ENSG00000091542 ALKBH5
load(file=r"(C:\Users\zxh\Desktop\R\meta collect\tcga_pcawg_rnaseq.Rdata)")
load(r"(C:\Users\zxh\Desktop\R\meta collect\ensg_symbol_relation.Rdata)")
TCGA_tpm_ensg_collapse=TCGA_tpm_ensg
TCGA_tpm_ensg_collapse[["KIDNEY"]]=do.call(cbind,list(TCGA_tpm_ensg_collapse[["KIRC"]],TCGA_tpm_ensg_collapse[["KICH"]],TCGA_tpm_ensg_collapse[["KIRP"]]))
TCGA_tpm_ensg_collapse[["KICH"]]=NULL
TCGA_tpm_ensg_collapse[["KIRP"]]=NULL
TCGA_tpm_ensg_collapse[["KIRC"]]=NULL
TCGA_tpm_ensg_collapse[["LAML"]]=NULL
TCGA_tpm_ensg_collapse[["DLBC"]]=NULL
TCGA_tpm_ensg_collapse[["COLORECTAAL"]]=do.call(cbind,list(TCGA_tpm_ensg_collapse[["COAD"]],TCGA_tpm_ensg_collapse[["PRAD"]]))
TCGA_tpm_ensg_collapse[["COAD"]]=NULL
TCGA_tpm_ensg_collapse[["PRAD"]]=NULL
library(doParallel) #
cl <- makeCluster(2)
registerDoParallel(cl)
foreach(i=1:33,.packages = c(library(survival),library(survminer),library(ggtext),library(cowplot))) %dopar% after_multiple_surviaval(genes="ENSG00000038427",rna.expr=TCGA_tpm_ensg_collapse,metaMatrix=pcawg_clinical_CDR,i) # returns list
stopCluster(cl)


lapply(choose_gene[[1]],function(x) multiple_surviaval(x))
lapply(choose_gene[[1]],function(x,dat) as.character(x))

if(F){
  lapply(choose_gene[[9]],function(x) multiple_surviaval(x))
  lapply(choose_gene[[9]],function(x,dat) as.character(x))
  
  library(doParallel) #
  cl <- makeCluster(1)
  registerDoParallel(cl)
  foreach(x=choose_gene[[9]],.packages = c(library(survival),library(survminer),library(ggtext),library(cowplot))) %dopar% multiple_surviaval(x) # returns list
  
  
  x="ENSG00000167306" 
}


for(i in 1:length(fusionALL1)){
  write.csv(fusionALL1[[i]], file=out_filePath[i],quote = T) 
}



outPath <- "./fusion gene lists/" ##输出路径
out_fileName <- sapply(names(fusionALL1),function(x){
  paste(x, ".csv", sep='')}) ##csv格式
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')}) ##输出路径名
##输出文件
for(i in 1:length(fusionALL1)){
  write.csv(fusionALL1[[i]], file=out_filePath[i],quote = T) 
}



outPath <- "./multiple survival/" ##输出路径
out_fileName <- sapply(names(fusionALL1),function(x){
  paste(x, ".pdf", sep='')}) ##csv格式
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')}) ##输出路径名
##输出文件


model <- coxph(s, data =dat)
summary(model,data=dat)
options(scipen=1)
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
with(new_dat,rcorr.cens(fp,Surv(time, event)  ))
# http://dni-institute.in/blogs/cox-regression-interpret-result-and-predict/
# ?surv_categorize
library(cowplot)
library(pheatmap)
if(T){
  library(survival)
  library(survminer)
  risk=dat
  risk$fp=fp
  res.cut <- surv_cutpoint(risk, time = "time", event = "event",
                           variables =c("fp"))
  summary(res.cut)
  plot(res.cut, "fp", palette = "npg")
  res.cat <- surv_categorize(res.cut)
  head(res.cat)
  library("survival")
  fit <- survfit(Surv(time, event) ~fp, data = res.cat)
  ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE)
  
}

risk=new_dat
risk$fp=fp
res.cut.va <- surv_cutpoint(risk, time = "time", event = "event",
                            variables =c("fp"))
summary(res.cut.va)
# res.cut.va$cutpoint$cutpoint=0.82#summary(res.cut.va)$cutpoint
res.cat.va <- surv_categorize(res.cut.va)
fit <- survfit(Surv(time, event) ~fp, data = res.cat.va)
ggsurvplot(fit, data = res.cat.va, risk.table = TRUE, conf.int = TRUE, pval=TRUE)
dev.off()
pdf("k1.pdf",height = 8,width = 8,onefile = F)
ggsurvplot(fit, data = res.cat.va, risk.table = TRUE, conf.int = TRUE, pval=TRUE)
dev.off()
library(cowplot)
save_plot()

# 按照一个cut成两个组
# library(Greg)
# test_data=subset(dat,select = c("OS","event","time","sample"))
# test_data=timeSplitter(test_data, 1825, 
#              time_var = "time",
#              # time_related_vars = c("age", "date"),
#              event_var = "event")
# test_data$event
# https://cran.r-project.org/web/packages/Greg/vignettes/timeSplitter.html