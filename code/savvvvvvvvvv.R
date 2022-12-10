#第二步，没有把v0转换为v0d，先看一看，看完全选一个，这个不作为画图步骤

load("kaiqiaole_v2.Rdata")
mat.out10$set="10"
mat.out5$set="5"
mat.out3$set="3"
mat.out3$number=sapply(stringr::str_split(mat.out3$X8,"\\+",simplify = F),function(x) length(x))
mat.out5$number=sapply(stringr::str_split(mat.out5$X8,"\\+",simplify = F),function(x) length(x))
mat.out10$number=sapply(stringr::str_split(mat.out10$X8,"\\+",simplify = F),function(x) length(x))
mat.out=unique(rbind(mat.out10,mat.out5,mat.out3))

# mat.out=subset(mat.out,number<50)
new_dat=v2
new_dat$event=as.numeric(new_dat$event)
new_dat$time=as.numeric(new_dat$time)
vvalidation=new_dat
vvalidation=subset(vvalidation,!is.na(vvalidation$margin0))
options(scipen=200)
library(Hmisc)
library(parallel)
cl <- makeCluster(getOption( "cl.cores" , 32))
clusterEvalQ(cl,c(library(survival),library(timeROC)))
clusterExport(cl, c("mat.out","v0","vvalidation"))#不加这个不认
mat.out[mat.out$number=="7",]
library(survival)
cc=parLapply(cl,c(lapply(mat.out$X8, function(x) as.formula(paste("Surv(time, event) ~",x)))),
               function(x)  1-with(vvalidation,Hmisc::rcorr.cens(predict(coxph(x, data =v0),vvalidation),survival::Surv(time, event)))["C Index"])
mat.out$X1=unlist(cc) 
mat.outv2=mat.out
# load("kaiqiaole_v1.Rdata")
# mat.outv1=mat.out

# kk=as.data.frame(c(mat.outv2$X1,mat.outv1$X2))

# fp <- predict(coxph(x, data =v0),vvalidation);
# options(scipen=200)
# cc4=1-with(new_dat,rcorr.cens(fp,survival::Surv(time, event)))["C Index"]


time_roc_resC[[1]] <- parLapply(cl,c(lapply(mat.out$X8, function(x) as.formula(paste("Surv(time, event) ~",x)))),
                           function(x) timeROC::timeROC(
                             T = v0$time,
                             delta = v0$event,
                             marker = predict(coxph(x, data =v0),v0),
                             cause = 1,
                             weighting="marginal",#"marginal",用的是km
                             times = c(365,730,1095,1460,1825 ),
                             ROC = TRUE,
                             iid = TRUE
                           )
)
time_roc_resC[[2]] <- parLapply(cl,c(lapply(mat.out$X8, function(x) as.formula(paste("Surv(time, event) ~",x,"+dcN01+margin0")))),
                                function(x) timeROC::timeROC(
                                  T = v0$time,
                                  delta = v0$event,
                                  marker = predict(coxph(x, data =v0),v0),
                                  cause = 1,
                                  weighting="marginal",#"marginal",用的是km
                                  times = c(365,730,1095,1460,1825 ),
                                  ROC = TRUE,
                                  iid = TRUE
                                )
)


time_roc_resCCC=lapply(1:2,function(i) data.frame(t(sapply(1:31,function(x) time_roc_resC[[i]][[x]]$AUC))))
time_roc_resCCC1=cbind(time_roc_resCCC[[1]],time_roc_resCCC[[2]])




time_roc_resV[[1]] <- parLapply(cl,c(lapply(mat.out$X8, function(x) as.formula(paste("Surv(time, event) ~",x)))),
                                function(x) timeROC::timeROC(
                                  T = vvalidation$time,
                                  delta = vvalidation$event,
                                  marker = predict(coxph(x, data =v0),vvalidation),
                                  cause = 1,
                                  weighting="marginal",#"marginal",用的是km
                                  times = c(365,730,1095,1460,1825 ),
                                  ROC = TRUE,
                                  iid = TRUE
                                )
)



time_roc_resV[[2]] <- parLapply(cl,c(lapply(mat.out$X8, function(x) as.formula(paste("Surv(time, event) ~",x,"+margin0")))),
                                function(x) timeROC::timeROC(
                                  T = vvalidation$time,
                                  delta = vvalidation$event,
                                  marker = predict(coxph(x, data =v0),vvalidation),
                                  cause = 1,
                                  weighting="marginal",#"marginal",用的是km
                                  times = c(365,730,1095,1460,1825 ),
                                  ROC = TRUE,
                                  iid = TRUE
                                )
)
time_roc_resV[[3]] <- parLapply(cl,c(lapply(mat.out$X8, function(x) as.formula(paste("Surv(time, event) ~",x,"+dcN01")))),
                                function(x) timeROC::timeROC(
                                  T = vvalidation$time,
                                  delta = vvalidation$event,
                                  marker = predict(coxph(x, data =v0),vvalidation),
                                  cause = 1,
                                  weighting="marginal",#"marginal",用的是km
                                  times = c(365,730,1095,1460,1825 ),
                                  ROC = TRUE,
                                  iid = TRUE
                                )
)
time_roc_resV[[4]] <- parLapply(cl,c(lapply(mat.out$X8, function(x) as.formula(paste("Surv(time, event) ~",x,"+dcN01+margin0")))),
                                function(x) timeROC::timeROC(
                                  T = vvalidation$time,
                                  delta = vvalidation$event,
                                  marker = predict(coxph(x, data =v0),vvalidation),
                                  cause = 1,
                                  weighting="marginal",#"marginal",用的是km
                                  times = c(365,730,1095,1460,1825 ),
                                  ROC = TRUE,
                                  iid = TRUE
                                )
)

time_roc_resVVV=lapply(1:4,function(i) data.frame(t(sapply(1:31,function(x) time_roc_resV[[i]][[x]]$AUC))))
time_roc_resVVV1=cbind(time_roc_resVVV[[1]],time_roc_resVVV[[2]],time_roc_resVVV[[3]],time_roc_resVVV[[4]])

mat.outv2zy=cbind(time_roc_resVVV1, mat.outv2,time_roc_resCCC1)            
# mat.outv2zy=cbind(time_roc_resVVV1, mat.outv2,data.frame(t(sapply(1:31,function(x) time_roc_resC[[x]]$AUC))))

write.csv(mat.outv2zy,"mat.outv2zy.csv")


time_roc_res <- parLapply(cl,c(lapply(mat.out$X8, function(x) as.formula(paste("Surv(time, event) ~",x))),   as.formula("Surv(time, event) ~fhpT01"), 
                          as.formula("Surv(time, event) ~dcN01"),
                          as.formula("Surv(time, event) ~stage02a"),
                          as.formula("Surv(time, event) ~margin0"),
                          as.formula("Surv(time, event) ~margin0+stage02a"),
                          as.formula("Surv(time, event) ~margin0+dcN01"),
                          as.formula("Surv(time, event) ~ENSG00000116285+ENSG00000170921+ENSG00000136870+ENSG00000103150+ENSG00000177042+ENSG00000187840+ENSG00000005884+dcN01")),
                          function(x) timeROC::timeROC(
                            T = vvalidation$time,
                            delta = vvalidation$event,
                            marker = predict(coxph(x, data =v0),vvalidation),
                            cause = 1,
                            weighting="marginal",#"marginal",用的是km
                            times = c(365,730,1095,1460,1825 ),
                            ROC = TRUE,
                            iid = TRUE
                          )
)
ENSG00000116285+ENSG00000170921+ENSG00000136870+ENSG00000103150+ENSG00000177042+ENSG00000187840+ENSG00000005884
ENSG00000116285+ENSG00000134369+ENSG00000170921+ENSG00000136870+ENSG00000103150+ENSG00000142188+ENSG00000177042+ENSG00000187840+ENSG00000005884+ENSG00000135245
ENSG00000116285+ENSG00000134369+ENSG00000170921+ENSG00000136870+ENSG00000103150+ENSG00000142188+ENSG00000198707+ENSG00000177042+ENSG00000187840+ENSG00000005884+ENSG00000148926+ENSG00000135245+ENSG00000164465
mat.out$number
roc.test=timeROC::compare(time_roc_res[[37]],time_roc_res[[38]]) #必须完成is.na()
roc.test
time_ROC_df <- lapply(1:length(time_roc_res),function(x) data.frame(time_roc_res[[x]]$TP[, 3],time_roc_res[[x]]$FP[, 3]))
stopCluster(cl)

legend=c("Risk score: high vs low","pT4 vs pT1-3","pN+ vs pN0","IIb-IV vs.I-IIa","R1-2 vs R0","Resection margin + Pathological stage","Resection margin + Pathological nodal stage")
legend=c(mat.out$number,"pT4 vs pT1-3","pN+ vs pN0","IIb-IV vs.I-IIa","R1-2 vs R0","Resection margin + Pathological stage","Resection margin + Pathological nodal stage","cb")
library(ggplot2)

pdf(file='time_roc_gene_compare_clinical_validation.pdf',height = 8,width = 8)
ggplot() +
  eval(substitute( lapply(1:length(time_roc_res), function(x) geom_line(aes(x =time_roc_res..x...FP...3.,  y = time_roc_res..x...TP...3. ),data=time_ROC_df[[x]], size = 0.8, color = x))))+
  eval(substitute( lapply(1:length(time_roc_res), function(x) annotate("text",x = 1, y = 0.71-0.1*x, size = 4,  label = paste0(legend[[x]]," ", sprintf("%.3f", time_roc_res[[x]]$AUC[[3]])), color = x, hjust = 1))))+
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 3) +
  theme_bw() +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  # scale_y_continuous(limits = c(0, 1))+
  #labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "plain", size = 11, color = "black"),
    axis.title.x = element_text(face = "plain", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "plain", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dev.off()
###########
v0$fp <- predict(coxph(as.formula(paste("Surv(time, event) ~","ENSG00000116285+ENSG00000170921+ENSG00000136870+ENSG00000103150+ENSG00000177042+ENSG00000187840+ENSG00000005884")), data =v0),v0)
res.cut <- surv_cutpoint(v0, time = "time", event = "event",variables =c("fp"))
summary(res.cut)
res.cut$cutpoint$cutpoint=-0.138
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(time, event) ~fp, data = res.cat,type  = "kaplan-meier",conf.type = "log-log")
summary(fit, times =c(365,730,1095,1460,1825))
ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, break.time.by = 365)

fit <- survfit(Surv(time, event) ~1, data = res.cat,type  = "kaplan-meier",conf.type = "log-log")
summary(fit, times =c(365,730,1095,1460,1825))
fit
####
vvalidation$fp <- predict(coxph(as.formula(paste("Surv(time, event) ~","ENSG00000116285+ENSG00000170921+ENSG00000136870+ENSG00000103150+ENSG00000177042+ENSG00000187840+ENSG00000005884")), data =v0),vvalidation)
res.cut <- surv_cutpoint(vvalidation, time = "time", event = "event",variables =c("fp"))
summary(res.cut)
res.cut$cutpoint$cutpoint=-0.138
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(time, event) ~fp, data = res.cat,type  = "kaplan-meier",conf.type = "log-log") #这么写等于sas
summary(fit, times =c(365,730,1095,1460,1825))
ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, break.time.by = 365)


fit <- survfit(Surv(time, event) ~1, data = res.cat,type  = "kaplan-meier",conf.type = "log-log")  #这么写等于sas
summary(fit, times =c(365,730,1095,1460,1825))
fit
