
dput(names(JDW_replustime))
JDW_replustime1111111=subset(JDW_replustime,select=c("Sex", "Age", "TS", "TSN", "size", "Lauren", "Borrmann", "Diff", "lyvF", "NF", "signet", "pT", "pN", 
  "pM", "eber", "Tumormargin", "surgerymargin", 
  "Surgicalapproach", "Surgicalway", "MOR", "NAC", "NAR", 
  "Smoking", "Drinking", "HD", "hypertension", "diabetes", "Hyperlipidemia", 
  "Antihypertensive", "antidiabetic", "Statins","sqts", "shts", "PAC", "PAR", "height", "weight", "BMI",  "neoadjuvant", "yqts", "dsts", "After2020", "shts3", "xjyqts3","xjdsts3"))
JDW_replustime1111111$BMI
fit_p <- glm(shts3~Age2,family=binomial(),data=box_completecase)
summary(fit_p)
wwwwwwwwww=summary(fit_p)
fit_p <- MASS::glm.nb(xjdsts~GGT+AST+TT+MCV+MPV,data=box)
s245::qdredge(fit_p,na.action = "na.omit")
s245::myQAIC(fit_p)
ww=step(fit_p,)

library(bbmle)

library(data.table)
library(Publish)
library(survival)
library(tidyr)
library(dplyr)
box$neoadjuvant=as.factor(box$neoadjuvant)
aaa=summary(regressionTable(fit_p))
aaa=as.data.frame(aaa)
openxlsx::write.xlsx(aaa,file ="aaa.xlsx")
plot(fit_p)
sub_pois <- subgroupAnalysis(fit_p,box,treatment="neoadjuvant",
                             subgroups=+`MCH`)

summary(regressionTable(fit_p))

# +offset(log(observationTime)

#########################################

library(doParallel) #并行处理包
cl <- makeCluster(30)#makeCluster(detectCores())
registerDoParallel(cl)

library(lars) 
library(glmnet) 
library(survival)
library(glmnet)
library(survival)
library(survminer)

lapply(as.data.frame(x), function(x) table(is.na(x)))
lapply(as.data.frame(x), function(x) is.numeric(x))
class(x)
# box$fhpT=ifelse(is.na(box$pT),NA,
#                 ifelse(box$pT%in%
box=JDW_replustime
box_completecase=box[complete.cases(box[,c("neoadjuvant","Sex", "Age", "TS", "size", "Lauren", "Borrmann", "Diff", "signet","pT","pN","pM","Smoking", "Drinking", "HD", "hypertension", "diabetes")]),]#,output_huayanzhibiao
box_completecase$TS=as.factor(box_completecase$TS)
box_completecase$Lauren=factor(box_completecase$Lauren,ordered = F)
box_completecase$Borrmann=factor(box_completecase$Borrmann,ordered = F)
box_completecase$Diff=factor(box_completecase$Diff,ordered = T)
box_completecase$signet=factor(box_completecase$signet,ordered = T)
box_completecase$Age2 = cut(box_completecase$Age,breaks = c(-Inf,65,Inf),right =T,ordered_result = T,labels = c(0,1))
box_completecase$pT1=ifelse(is.na(box_completecase$pT),NA,
                            ifelse(box_completecase$pT=="0","0",
                                   ifelse(box_completecase$pT%in%c("1","1a","1b"),"1",
                                          ifelse(box_completecase$pT%in%c("2","2a","2b"),"2",
                                                 ifelse(box_completecase$pT%in%c("3","3a","3b"),"3",
                                                        ifelse(box_completecase$pT%in%c("4","4a"),"4",
                                                               ifelse(box_completecase$pT=="4b","5",NA)))))))
table(box_completecase$pT1)
box_completecase$pN1=ifelse(is.na(box_completecase$pN),NA,
                            ifelse(box_completecase$pN=="0","0",
                                   ifelse(box_completecase$pN%in%c("1","1a","1b"),"1",
                                          ifelse(box_completecase$pN%in%c("2","2a","2b"),"2",
                                                 ifelse(box_completecase$pN%in%c("3","3a"),"3",
                                                        ifelse(box_completecase$pN =="3b","5",NA))))))
box_completecase$pN1_correct=ifelse(is.na(box_completecase$pN),NA,
                            ifelse(box_completecase$pN=="0","0",
                                   ifelse(box_completecase$pN%in%c("1","1a","1b"),"1",
                                          ifelse(box_completecase$pN%in%c("2","2a","2b"),"2",
                                                 ifelse(box_completecase$pN%in%c("3","3a"),"3",
                                                        ifelse(box_completecase$pN =="3b","4",NA))))))
box_completecase$tplusn=as.character(as.numeric(box_completecase$pT1)+as.numeric(box_completecase$pN1))
# box_completecase$pstage3 = ifelse(box_completecase$tplusn==6 & box_completecase$pT==4 , "IIIA",NA)
box_completecase$pstage = ifelse(
  box_completecase$pM == "1","8",
  ifelse(box_completecase$tplusn =="1", "1",
         ifelse(box_completecase$tplusn =="2", "2",
                ifelse(box_completecase$tplusn =="3", "3",
                       ifelse(box_completecase$tplusn =="4", "4",
                              ifelse(box_completecase$tplusn =="5", "5",
                                     ifelse(box_completecase$tplusn=="6" & box_completecase$pT1=="4" , "5",
                                            ifelse(box_completecase$tplusn%in%c("6","7"), "6", 
                                                   ifelse(box_completecase$tplusn%in%c("8","9","10"), "7", NA)))))))))
box_completecase$pstage1=ifelse(is.na(box_completecase$pstage),NA,
                                ifelse(box_completecase$pstage %in%c("1", "2"),"1",
                                       ifelse(box_completecase$pstage %in%c("3","4"), "2",
                                              ifelse(box_completecase$pstage %in%c("5","6","7"), "3",
                                                     ifelse(box_completecase$pstage =="8", "4",NA)))))
box_completecase$pT12=ifelse(is.na(box_completecase$pT),NA,
                           ifelse(box_completecase$pT%in%c("0","1","1a","1b"),"0",
                                  ifelse(box_completecase$pT%in%c("2","2a","2b",
                                                                  "3","3a","3b","4","4a","4b"),"1",NA)))
box_completecase$pN01=ifelse(is.na(box_completecase$pN),NA,
                           ifelse(box_completecase$pN=="0","0",
                                  ifelse(box_completecase$pN%in%c("1","1a","1b","2",
                                                                  "2a","2b","3","3a","3b"),"1",NA)))
box_completecase$pM=as.numeric(box_completecase$pM)
dput(names(box_completecase))
box_completecase$Drinking
y=box_completecase$shts3
box_completecase= box_completecase %>%mutate(across(where(is.character),as.numeric)) 

x=box_completecase[,c("neoadjuvant", "Sex", "Age", "TS", "size", "Lauren", "Borrmann", "Diff", "signet", "pM", "Smoking", "Drinking", "HD", "hypertension", "diabetes", "Age2", "pT1", "pN1_correct","pstage", "pstage1", "pT12", "pN01")]#,output_huayanzhibiao
x=Hmisc::asNumericMatrix(x)
nrow(na.omit(x))
box_completecase_surgeryonly=subset(box_completecase,box_completecase$neoadjuvant==1)
box_completecase_neoadjuvant=subset(box_completecase,box_completecase$neoadjuvant==0)
write.csv(x,file="bbbb.csv")
lapply(as.data.frame(x),function(x) class(x))

FGT3 <- function(i){
  library(survival)
  set.seed(i)
  cv_fit <- glmnet::cv.glmnet(x, y, alpha=1,type.measure="deviance",standardize = TRUE,family = 'binomial',nfolds =3)
  fit <- glmnet::glmnet(x=x, y=y, alpha = 1,family ='binomial',standardize = TRUE,lambda=cv_fit$lambda.min)#lambda.1se
  choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
  if(length(choose_gene)!=0){s=as.formula(paste("y ~" ,paste0("`",choose_gene,"`",collapse = "+")))
  model <- glm(formula=s,family = binomial(), data =box_completecase)
  #
  new_dat=box_completecase
  fp <- predict(model,new_dat) ;
  library(Hmisc)
  options(scipen=200)
  cc1=with(new_dat,rcorr.cens(fp,y))["C Index"]
  if(cc1<0.5){
  cc1=1-with(new_dat,rcorr.cens(fp,y))["C Index"]
  }
  new_dat=box_completecase_surgeryonly
  fp <- predict(model,new_dat) ;
  library(Hmisc)
  options(scipen=200)
  cc2=with(new_dat,rcorr.cens(fp,new_dat$shts3))["C Index"]
  if(cc2<0.5){
    cc2=1-with(new_dat,rcorr.cens(fp,new_dat$shts3))["C Index"]
  }
  new_dat=box_completecase_neoadjuvant
  fp <- predict(model,new_dat) ;
  library(Hmisc)
  options(scipen=200)
  cc3=with(new_dat,rcorr.cens(fp,new_dat$shts3))["C Index"]
  if(cc3<0.5){
    cc3=1-with(new_dat,rcorr.cens(fp,new_dat$shts3))["C Index"]
  }
  mat.out=list(c(cc1,cc2,cc3,paste(choose_gene,collapse = "+")))}else{ mat.out=list(rep(NA,3+1))}
  mat.out
}
mat=foreach(i=1:1000,.combine = "c")%dopar% FGT3(i)
mat.out3 <- data.frame(matrix(unlist(mat), nrow=length(mat), byrow=T),stringsAsFactors=FALSE)
FGT5 <- function(i){
  library(survival)
  set.seed(i)
  cv_fit <- glmnet::cv.glmnet(x, y, alpha=1,type.measure="deviance",standardize = TRUE,family = 'binomial',nfolds =5)
  fit <- glmnet::glmnet(x=x, y=y, alpha = 1,family ='binomial',standardize = TRUE,lambda=cv_fit$lambda.min)#lambda.1se
  choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
  if(length(choose_gene)!=0){s=as.formula(paste("y ~" ,paste0("`",choose_gene,"`",collapse = "+")))
  model <- glm(formula=s,family = binomial(), data =box_completecase)
  #
  new_dat=box_completecase
  fp <- predict(model,new_dat) ;
  library(Hmisc)
  options(scipen=200)
  cc1=with(new_dat,rcorr.cens(fp,y))["C Index"]
  if(cc1<0.5){
    cc1=1-with(new_dat,rcorr.cens(fp,y))["C Index"]
  }
  new_dat=box_completecase_surgeryonly
  fp <- predict(model,new_dat) ;
  library(Hmisc)
  options(scipen=200)
  cc2=with(new_dat,rcorr.cens(fp,new_dat$dsts3))["C Index"]
  if(cc2<0.5){
    cc2=1-with(new_dat,rcorr.cens(fp,new_dat$dsts3))["C Index"]
  }
  new_dat=box_completecase_neoadjuvant
  fp <- predict(model,new_dat) ;
  library(Hmisc)
  options(scipen=200)
  cc3=with(new_dat,rcorr.cens(fp,new_dat$dsts3))["C Index"]
  if(cc3<0.5){
    cc3=1-with(new_dat,rcorr.cens(fp,new_dat$dsts3))["C Index"]
  }
  mat.out=list(c(cc1,cc2,cc3,paste(choose_gene,collapse = "+")))}else{ mat.out=list(rep(NA,3+1))}
  mat.out
}
mat=foreach(i=1:1000,.combine = "c")  %dopar% FGT5(i)
mat.out5 <- data.frame(matrix(unlist(mat), nrow=length(mat), byrow=T),stringsAsFactors=FALSE)
FGT10 <- function(i){
  library(survival)
  set.seed(i)
  cv_fit <- glmnet::cv.glmnet(x, y, alpha=1,type.measure="deviance",standardize = TRUE,family = 'binomial',nfolds =10)
  fit <- glmnet::glmnet(x=x, y=y, alpha = 1,family ='binomial',standardize = TRUE,lambda=cv_fit$lambda.min)#lambda.1se
  choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
  if(length(choose_gene)!=0){s=as.formula(paste("y ~" ,paste0("`",choose_gene,"`",collapse = "+")))
  model <- glm(formula=s,family = binomial(), data =box_completecase)
  #
  new_dat=box_completecase
  fp <- predict(model,new_dat) ;
  library(Hmisc)
  options(scipen=200)
  cc1=with(new_dat,rcorr.cens(fp,y))["C Index"]
  if(cc1<0.5){
    cc1=1-with(new_dat,rcorr.cens(fp,y))["C Index"]
  }
  new_dat=box_completecase_surgeryonly
  fp <- predict(model,new_dat) ;
  library(Hmisc)
  options(scipen=200)
  cc2=with(new_dat,rcorr.cens(fp,new_dat$dsts3))["C Index"]
  if(cc2<0.5){
    cc2=1-with(new_dat,rcorr.cens(fp,new_dat$dsts3))["C Index"]
  }
  new_dat=box_completecase_neoadjuvant
  fp <- predict(model,new_dat) ;
  library(Hmisc)
  options(scipen=200)
  cc3=with(new_dat,rcorr.cens(fp,new_dat$dsts3))["C Index"]
  if(cc3<0.5){
    cc3=1-with(new_dat,rcorr.cens(fp,new_dat$dsts3))["C Index"]
  }
  mat.out=list(c(cc1,cc2,cc3,paste(choose_gene,collapse = "+")))}else{ mat.out=list(rep(NA,3+1))}
  mat.out
}
mat=foreach(i=1:1000,.combine = "c")  %dopar% FGT10(i)
mat.out10 <- data.frame(matrix(unlist(mat), nrow=length(mat), byrow=T),stringsAsFactors=FALSE)
####3###
mat.out=unique(rbind(mat.out10,mat.out5,mat.out3))
mat.out$number=sapply(stringr::str_split(mat.out$X4,"\\+",simplify = F),function(x) length(x))
stopCluster(cl)


library(lars) 
library(glmnet) 
library(survival)
library(glmnet)
#check
set.seed(4)
cv_fit <- glmnet::cv.glmnet(x, y, alpha=1,type.measure="deviance",standardize = TRUE,family = 'binomial',nfolds =3)
pdf(file =paste0("Figure 4A.dev.pdf"),height = 10,width = 10)
plot(cv_fit)
dev.off()
pdf(file =paste0("Figure 4B.rlambda.pdf"),height = 10,width = 10)
set.seed(4)
model_lasso <- glmnet(x, y, alpha=1,family ='binomial', nlambda=10000)
print(model_lasso)
# 列%Dev代表了由模型解释的残差的比例，对于线性模型来说就是模型拟合的R^2(R-squred)。
# 它在0和1之间，越接近1说明模型的表现越好，
# 如果是0，说明模型的预测结果还不如直接把因变量的均值作为预测值来的有效。
# 使用area under the ROC curve, CV 选择压缩参数lambda
# 再设置一次set.seed
plotmo::plot_glmnet(model_lasso,s=cv_fit$lambda.min, xvar = "rlambda", label = TRUE)
dev.off()
set.seed(4)
fit <- glmnet::glmnet(x=x, y=y, alpha = 1,family ='binomial',standardize = TRUE,lambda=cv_fit$lambda.min)#lambda.1se
print(fit)
# fit the model
# head(coef(model_lasso, s=c(model_lasso$lambda[282],0.24690)))#??????
# 两条虚线分别指示了两个特殊的λ值:
c(log(cv_fit$lambda.min),log(cv_fit$lambda.1se))
c(cv_fit$lambda.min,cv_fit$lambda.1se)

choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
choose_gene=names(coeff$ipf)[-1]
s=as.formula(paste("y ~" ,paste0("`",choose_gene,"`",collapse = "+")))
model <- glm(formula=s,family = binomial(), data =box_completecase)
#
new_dat=box_completecase
fp <- predict(model,new_dat) ;
library(Hmisc)
options(scipen=200)
cc1=with(new_dat,rcorr.cens(fp,y))["C Index"]
if(cc1<0.5){
  cc1=1-with(new_dat,rcorr.cens(fp,y))["C Index"]
}
cc1
# k <- dim(x)[1]
# predictions <- c()
# for (i in 1:k) {
#   model <- glmnet(x[-i,], y[-i], family="binomial")
#   predictions <- c(predictions, predict(model, newx=x[i,]))
# }
# library(pROC)
# roc(y, predictions)

library(data.table)
get_coe <- function(the_fit,the_lamb){
  Coefficients <- coef(the_fit, s = the_lamb)
  Active.Index <- which(Coefficients != 0)
  Active.Coefficients <- Coefficients[Active.Index]
  re <- data.frame(rownames(Coefficients)[Active.Index],Active.Coefficients)
  re <- data.table('var_names'=rownames(Coefficients)[Active.Index],
                   'coef'=Active.Coefficients)
  # 计算系数的指数次方，表示x每变动一个单位对y的影响倍数
  re$expcoef <- exp(re$coef)
  return(re[order(expcoef)])
}
get_coe(cv_fit,cv_fit$lambda.min)

get_plot<- function(the_fit,the_fit_cv,the_lamb,toplot = seq(1,50,2)){
  Coefficients <- coef(the_fit, s = the_lamb)
  Active.Index <- which(Coefficients != 0)
  coeall <- coef(the_fit, s = the_fit_cv$lambda[toplot])
  coe <- coeall[Active.Index[-1],]
  ylims=c(-max(abs(coe)),max(abs(coe)))
  sp <- spline(log(the_fit_cv$lambda[toplot]),coe[1,],n=100)
  plot(sp,type='l',col =1,lty=1,
       ylim = ylims,ylab = 'Coefficient', xlab = 'log(lambda)')
  abline(h=0)
  for(i in c(2:nrow(coe))){
    lines(spline(log(the_fit_cv$lambda[toplot]),coe[i,],n=1000),
          col =i,lty=i)
  }
  legend("bottomright",legend=rownames(coe),col=c(1:nrow(coe)),
         lty=c(1:nrow(coe)),
         cex=0.5)
}
# 传入最优lambda-1，从而保留更多变量
get_plot(model_lasso,cv_fit,exp(log(cv_fit$lambda.min)-1))

#这些只能用来logistic
library(ROCR)
get_confusion_stat <- function(pred,y_real,threshold=0.5){
  # auc
  tmp <- prediction(as.vector(pred),y_real)
  auc <- unlist(slot(performance(tmp,'auc'),'y.values'))
  # statistic
  pred_new <- as.integer(pred>threshold)
  tab <- table(pred_new,y_real)
  if(nrow(tab)==1){
    print('preds all zero !')
    return(0)
  }
  TP <- tab[2,2]
  TN <- tab[1,1]
  FP <- tab[2,1]
  FN <- tab[1,2]
  accuracy <- round((TP+TN)/(TP+FN+FP+TN),4)
  recall_sensitivity <- round(TP/(TP+FN),4)
  precision <- round(TP/(TP+FP),4)
  specificity <- round(TN/(TN+FP),4)
  # 添加，预测的负例占比（业务解释：去除多少的样本，达到多少的recall）
  neg_rate <- round((TN+FN)/(TP+TN+FP+FN),4)
  re <- list('AUC' = auc,
             'Confusion_Matrix'=tab,
             'Statistics'=data.frame(value=c('accuracy'=accuracy,
                                             'recall_sensitivity'=recall_sensitivity,
                                             'precision'=precision,
                                             'specificity'=specificity,
                                             'neg_rate'=neg_rate)))
  return(re)
}
# 结合数据预处理部分，得到模型在测试集上的表现
get_eval <- function(expr,y_name,theta=0.5,the_fit=fit,the_lamb=fit_cv$lambda.min){
  y <- y_name
  # x <- t(expr)# zxh 注意啥时候需要改的 20220306
  x <- as.matrix(x)
  pred <- predict(the_fit,newx=x,s=the_lamb,type = 'response')#the type of prediction required. The default is on the scale of the linear predictors; the alternative "response" is on the scale of the response variable. Thus for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives the predicted probabilities. The "terms" option returns a matrix giving the fitted values of each term in the model formula on the linear predictor scale.
  print(get_confusion_stat(pred,y, theta))
}
get_eval(expr=x,
         y_name=y,
         the_fit=fit,
         the_lamb=cv_fit$lambda.min,
         )


#load libraries
library(data.table)
library(Publish)
library(survival)
library(tidyr)
library(dplyr)

model <- glm(formula=s,family = binomial(), data =box_completecase)




traceR=as.data.frame(box_completecase)
# traceR$dsts3=as.factor(traceR$dsts3)
dput(names(traceR))
traceR1=subset(traceR,select=c("病案号", "neoadjuvant", "Sex", "Age", "TS", "size", "Lauren", 
  "Borrmann", "Diff", "signet", "pT", "pN", "pM", "Smoking", "Drinking", 
  "HD", "hypertension", "diabetes", "Age2","pT1", "pN1", "pN1_correct", "tplusn", "pstage", "pstage1", "pT12", 
  "pN01"))
traceR1= traceR1 %>%mutate(across(where(is.numeric),as.factor)) 
traceR=cbind(traceR1,traceR[,c("dsts3",output_huayanzhibiao)])
traceR$dsts3=factor(traceR$dsts3)
fit_log <- glm(dsts3 ~ neoadjuvant + Age + Sex + ALT + AST + `CK-MB` + CRE + 
                 G + GGT + GLU + HBDH + PALB + TP + `EOS#` + `EOS%` + 
                 `MONO#` + `MONO%` + MPV + `P-LCR` + PDW + 
                 SOD + IBIL + TRANSFE + HCT + ALB + ALP + IGG + LDH + TBA + 
                 MCH + MCV + RBC + `RDW-CV` + `RDW-SD` + WBC + 
                 HB + UREA + DBIL + CK,family="binomial",data=traceR)
sub_log <- subgroupAnalysis(fit_log,traceR,treatment="dsts3",
                            subgroups=~Sex, factor.reference="inline")
traceR$Sex
uwlaalal=plot(sub_pois)
plot(fit_log)
x=sub_log

plot.subgroupAnalysis <- function(x,...)
{
  if (class(x)[1]!="subgroupAnalysis") stop("Object not of class subgroupAnalysis")
  num <- length(names(x))
  plotcols<-x[,(num-4):(num-2)]
  tabcols <-x[,c(1:2,num)]
  Publish::plotConfidence(x=plotcols, labels=tabcols)
}