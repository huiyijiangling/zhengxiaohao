# ### ---------------
# ###
# ### Create: Jianming Zeng
# ### Date: 2019-04-02 21:59:01
# ### Email: jmzeng1314@163.com
# ### Blog: http://www.bio-info-trainee.com/
# ### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
# ### CAFS/SUSTC/Eli Lilly/University of Macau
# ### Update Log:   2019-04-02  second version
# ###
# ### ---------------
# 
# # https://mp.weixin.qq.com/s/J-vaFq1Vv-zR1LC0Wop1zw
# # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5051943/figure/djw200-F2/
# # https://www.ncbi.nlm.nih.gov/pubmed/24893932
# # Sixteen of the 111 most significantly altered miRNAs were associated with OS across different clinical subclasses of the TCGA-derived STAD cohort. 
# # A linear prognostic model of eight miRNAs (miR-31, miR-196b, miR-766, miR-519a-1, miR-375, miR-187, miR-331 and miR-101-1) was constructed and weighted by the importance scores from the supervised principal component method to divide patients into high- and low-risk groups. 
# # Patients assigned to the high-risk group exhibited poor OS compared with patients in the low-risk group (hazard ratio [HR]=1.99, P <0.001). 
# # The eight-miRNA signature is an independent prognostic marker of OS of STAD patients and demonstrates good performance for predicting 5-year OS (Area Under the respective ROC Curves [AUC] = 0.626, P = 0.003), especially for non-smokers (AUC = 0.686, P = 0.023).
# 
# rm(list=ls())
# options(stringsAsFactors = F)
# 
# Rdata_dir='./Rdata/'
# Figure_dir='./figures/'
# 
# 
# library(survival)
# library(survminer)
# 
# # 这里举例的文章不一样，所以不再使用前面步骤的数据。
# 
# # 而是对TCGA-STAD-miRNA重新处理，正好跟前面的步骤对应学习。
# 
# # if(F){
# 
#   library(RTCGA.clinical) 
#   ??RTCGA.clinical
#   meta <- STAD.clinical
#   tmp=as.data.frame(colnames(meta))
#   meta[(grepl('patient.bcr_patient_barcode',colnames(meta)))]
#   meta[(grepl('patient.days_to_last_followup',colnames(meta)))]
#   meta[(grepl('patient.days_to_death',colnames(meta)))]
#   meta[(grepl('patient.vital_status',colnames(meta)))]
#   ## patient.race  # patient.age_at_initial_pathologic_diagnosis # patient.gender 
#   # patient.stage_event.clinical_stage
#   meta=as.data.frame(meta[c('patient.bcr_patient_barcode','patient.vital_status',
#                             'patient.days_to_death','patient.days_to_last_followup',
#                             'patient.race',
#                             'patient.age_at_initial_pathologic_diagnosis',
#                             'patient.gender' ,
#                             'patient.stage_event.pathologic_stage')])
#   #meta[(grepl('patient.stage_event.pathologic_stage',colnames(meta)))]
#   load(file="STAD_GDCRNATOOLS_ebv.Rdata") 
# expr=rnaCounts
# group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
# table(group_list)
# # 
# # if(F){
#   exprSet=na.omit(expr)
#   exprSet=exprSet[,group_list=='tumor']
#   
#   head(meta)
#   colnames(meta)
#   # write.csv(meta,file="meta.csv")
#   # write.csv(clinicalDa,file="STAD_clinical.csv")
#   meta=read.csv(file="meta.csv",stringsAsFactors = F,sep=",",fill = TRUE,encoding = "UTF-8",header=T)
# 
#   meta[,3][is.na(meta[,3])]=0
#   meta[,4][is.na(meta[,4])]=0
#   meta$days=as.numeric(meta[,3])+as.numeric(meta[,4])
#   meta=meta[,c(1:2,5:9)]
#   colnames(meta)
#   colnames(meta)=c('ID','event','race','age','gender','stage',"days") 
#   library(survival)
#   library(survminer)
#   meta$event=ifelse(meta$event=='Alive',0,1)
#   meta$age=as.numeric(meta$age)
#   library(stringr) 
#   meta$stage=str_split(meta$stage,' ',simplify = T)[,2]
#   table(  meta$stage)
#   boxplot(meta$age)
#   meta$age_group=ifelse(meta$age>median(meta$age,na.rm = T),'older','younger')
#   table(meta$race)
#   meta$time=meta$days/30
#   phe=meta
#   
#   head(phe)
#   phe$ID=toupper(phe$ID) 
#   exprSet=rnaCounts[,1:23]
#   phe=phe[match(substr(colnames(exprSet),1,12),phe$ID),]
#   head(phe)
#   exprSet[1:4,1:4]
# 
#   save(exprSet,phe,file=file.path(Rdata_dir,'TCGA-STAD-survival_input.Rdata'))
# #  }
# # 
#   mySurv=with(phe,Surv(time, event))
#   log_rank_p <- apply(exprSet , 1 , function(gene){
#     # gene=exprSet[1,]
#     phe$group=ifelse(gene>median(gene),'high','low')  
#     data.survdiff=survdiff(mySurv~group,data=phe)
#     p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
#     return(p.val)
#   })
# # load(file=file.path(Rdata_dir,'TCGA-STAD-survival_input.Rdata'))
# head(phe)
# 
# # 这个时候是515个病人的673个miRNA表达矩阵。
# 
# ## 挑选感兴趣的基因构建coxph模型 
# 
# # miR-31, miR-196b, miR-766, miR-519a-1, miR-375, miR-187, miR-331 and miR-101-1
# # hsa-mir-31,hsa-mir-196b,hsa-mir-766,hsa-mir-519a-1,hsa-mir-375,hsa-mir-187,hsa-mir-331,hsa-mir-101-1
# # e=t(exprSet[c('ENSG00000136108')])
# # e=log2(e+1)
# # colnames(e)=c('ENSG00000136108')
# # 从这里开始关闭
quantile(dat$age,probs = seq(0, 1, 0.05))
?quantile
exprSet3=as.data.frame(as.matrix(t(exprSet_quant_with_clinical)))

exprSet3$bah=rownames(exprSet3)
dat=merge(phe,exprSet3,by="bah")
# dat$gender=factor(dat$gender)
# dat$stage=factor(dat$stage)
library(survival)
library(survminer)
dat$event=as.numeric(dat$RFI)#dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$RFI.time)#dat$time=as.numeric(dat$time)
s=Surv(time, event) ~  age#NerveF+gulijiejie+hualiao+rbwz3
#KIF15+FEN1+TTF2+MSI2+KYNU+ZNF562+ACLY+CXCR4+KIF21B+SLC12A7+ZNF823+NHS  

model <- coxph(s, data =dat)
summary(model,data=dat)
options(scipen=1)

ggforest(model, data =dat, 
         main = "Forest plot (Hazard ratio)", 
         cpositions = c(0.10, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)+
ggsave('Forest.pdf',width = 8,height=6,dpi = 600)

# 对于生存数据预后模型评价很多采用C-index ，但c-index展示没有roc曲线图来的直观
new_dat=dat
new_dat$event=as.numeric(new_dat$event)
new_dat$time=as.numeric(new_dat$time)
# ??predict.coxph
# ??predict
# fp <- predict(model,new_dat,type="risk");boxplot(fp)#riskscore e^(lp的值)
# fp <- predict(model,new_dat,type="expected") ;boxplot(fp)#负的
# fp <- predict(model,new_dat,type="term") ;boxplot(fp)#负的
# fp <- predict(model,new_dat,type="survival") ;boxplot(fp)#负的
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
  # res.cut <- surv_cutpoint(dat, time = "time", event = "event",
  #                          variables =c('DCLRE1B','CBFB','STRIP2','MSI2','KYNU','DUSP14','PSMD12','ACLY','KIF21B','APCDD1','MT1H'))
  # summary(res.cut)
  risk=dat
  risk$fp=fp
  # intercept=-sum(coef(model)*apply(kkkkk,2,mean))
  # risk$fp=risk$fp-intercept
  res.cut <- surv_cutpoint(risk, time = "time", event = "event",
                           variables =c("fp"))
  summary(res.cut)
  # palette = "npg" (nature publishing group), see ?ggpubr::ggpar
  plot(res.cut, "fp", palette = "npg")
  ## $DEPDC1
  res.cat <- surv_categorize(res.cut)
  head(res.cat)
  library("survival")
  fit <- survfit(Surv(time, event) ~fp, data = res.cat)
  ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE)
  
}
#
risk=new_dat
risk$fp=fp
res.cut.va <- surv_cutpoint(risk, time = "time", event = "event",
                           variables =c("fp"))
summary(res.cut.va)
res.cut.va$cutpoint$cutpoint=0.82#summary(res.cut.va)$cutpoint
res.cat.va <- surv_categorize(res.cut.va)
fit <- survfit(Surv(time, event) ~fp, data = res.cat.va)
ggsurvplot(fit, data = res.cat.va, risk.table = TRUE, conf.int = TRUE, pval=TRUE)
dev.off()
pdf("k1.pdf",height = 8,width = 8,onefile = F)
ggsurvplot(fit, data = res.cat.va, risk.table = TRUE, conf.int = TRUE, pval=TRUE)
dev.off()
# https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
summary(res.cut.va)$cutpoint
fp_dat=data.frame(s=1:length(fp),
                  v=as.numeric(sort(fp )),
                  risk=ifelse(as.numeric(sort(fp ))>res.cut.va$cutpoint$cutpoint,1,0))
sur_dat=data.frame(s=1:length(fp),#sort(fp),#也可以
                   t=new_dat[names(sort(fp )),'time'] ,
                   Status=new_dat[names(sort(fp)),'event']  ) 
sur_dat$Status=ifelse(sur_dat$Status==0,'Alive','Death')
exp_dat=new_dat[names(sort(fp )),c('KIF15','FEN1','ZFP69B','SP6','SPARC','TTF2','MSI2','KYNU','ACLY','KIF21B','SLC12A7','ZNF823')]
plot.point=ggplot(fp_dat,aes(x=s,y=v))+geom_point(aes(col=risk))+ xlab(" ") + ylab("Risk");print(plot.point)
plot.sur=ggplot(sur_dat,aes(x=s,y=t))+geom_point(aes(col=Status))+theme(legend.background = element_blank())+ xlab(" ") + ylab("Time(months)");print(plot.sur)
#theme(legend.position=c(0.95, 0.8),legend.background = element_blank()) 这样能把图例放进去
mycolors <- colorRampPalette(c("green","black","red"), bias = 1.2)(100)
# tmp=t(scale(exp_dat))
# tmp[tmp > 1] = 1
# tmp[tmp < -1] = -1
tmp=t(exp_dat)
tmp[tmp > 2] = 2
tmp[tmp < -2] = -2
plot.h=pheatmap(tmp,show_colnames = F,cluster_cols = T,onefile = T)#
plot.h=pheatmap(main = " ",tmp,show_colnames = F,cluster_cols = F,col= mycolors)#,col= mycolors
plot_grid(plot.point, plot.sur, plot.h$gtable,
          labels = c("A", "B","C"),
          axis = "tb",
          align = 'v',ncol = 1)
dev.off()
pdf("Construction_risk_score.pdf",width = 6,height=10)
plot_grid(plot.point, plot.sur, plot.h$gtable,
          labels = c("A", "B","C"),
          axis = "tb",
          align = 'v',ncol = 1)
dev.off()
#exprSet3 
coef(model)
kkkkk=new_dat[,c("A1","A20","A36","A42","A54","A65","A98","A105","A106","A111","A114","A117")]
coef(model)*apply(kkkkk,2,median)
as.data.frame(coef(model))
sum(coef(model)*kkkkk[1,])
sum(coef(model)*apply(kkkkk,2,median))
sum(coef(model)*apply(kkkkk,2,mean))
intercept=-sum(coef(model)*apply(kkkkk,2,mean))
intercept
sum(coef(model)*kkkkk[1,])-sum(coef(model)*apply(kkkkk,2,mean))
risk$fp=risk$fp-intercept
