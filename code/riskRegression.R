# binary outcome
library(riskRegression)
library(lava)
set.seed(18)
learndat <- sampleData(48,outcome="binary")
testdat <- sampleData(40,outcome="binary")

## score logistic regression models
lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
lr2 = glm(Y~X3+X5,data=learndat,family=binomial)
Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5)"=lr2),formula=Y~1,data=testdat)

## ROC curve and calibration plot
xb=Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5+X6)"=lr2),formula=Y~1,
         data=testdat,plots=c("calibration","ROC"))
## Not run: plotROC(xb)
plotCalibration(xb)

## End(Not run)

## compute AUC for a list of continuous markers
markers = as.list(testdat[,.(X6,X7,X8,X9,X10)])
Score(markers,formula=Y~1,data=testdat,metrics=c("auc"))

# cross-validation
## Not run: 
learndat=sampleData(400,outcome="binary")
lr1a = glm(Y~X6,data=learndat,family=binomial)
lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
## bootstrap cross-validation
x1=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=100)
x1
## leave-one-out and leave-pair-out bootstrap
x2=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,
         split.method="loob",
         B=100,plots="calibration")
x2

## End(Not run)
# survival outcome

# Score Cox regression models
## Not run: library(survival)
library(rms)
library(prodlim)
set.seed(18)
trainSurv <- sampleData(100,outcome="survival")
testSurv <- sampleData(40,outcome="survival")
cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
         formula=Surv(time,event)~1,data=testSurv,conf.int=FALSE,times=c(5,8))
xs

## End(Not run)

# Integrated Brier score
## Not run: 
xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
         formula=Surv(time,event)~1,data=testSurv,conf.int=FALSE,
         summary="ibs",
         times=sort(unique(testSurv$time)))

## End(Not run)

# time-dependent AUC for list of markers
## Not run: survmarkers = as.list(testSurv[,.(X6,X7,X8,X9,X10)])
Score(survmarkers,
      formula=Surv(time,event)~1,metrics="auc",data=testSurv,
      conf.int=TRUE,times=c(5,8))

# compare models on test data
Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
      formula=Surv(time,event)~1,data=testSurv,conf.int=TRUE,times=c(5,8))

## End(Not run)
# crossvalidation models in traindata
## Not run: 
library(survival)
set.seed(18)
trainSurv <- sampleData(400,outcome="survival")
cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
x1 = Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
           formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=c(5,8),
           split.method="loob",B=100,plots="calibration")

x2= Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
          formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=c(5,8),
          split.method="bootcv",B=100)

## End(Not run)

# restrict number of comparisons
## Not run: 
Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
      formula=Surv(time,event)~1,data=trainSurv,contrasts=TRUE,
      null.model=FALSE,conf.int=TRUE,times=c(5,8),split.method="bootcv",B=3)

# competing risks outcome
set.seed(18)
trainCR <- sampleData(40,outcome="competing.risks")
testCR <- sampleData(40,outcome="competing.risks")
library(riskRegression)
library(cmprsk)
# Cause-specific Cox regression
csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR)
csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR)
# Fine-Gray regression
fgr1 = FGR(Hist(time,event)~X1+X2+X7+X9,data=trainCR,cause=1)
fgr2 = FGR(Hist(time,event)~X3+X5+X6,data=trainCR,cause=1)
Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2,
           "FGR(X1+X2+X7+X9)"=fgr1,"FGR(X3+X5+X6)"=fgr2),
      formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(5,8))

## End(Not run)



## Not run: 
# reproduce some results of Table IV of Blanche et al. Stat Med 2013
data(Paquid)
ResPaquid <- Score(list("DSST"=-Paquid$DSST,"MMSE"=-Paquid$MMSE),
                   formula=Hist(time,status)~1,
                   data=Paquid,
                   null.model = FALSE,
                   conf.int=TRUE,
                   metrics=c("auc"),
                   times=c(3,5,10),
                   plots="ROC")
ResPaquid
plotROC(ResPaquid,time=5)
Paquid$DSST[1]=1
## End(Not run)
## Not run: 
# parallel options
# by erikvona: Here is a generic example of using future
# and doFuture, works great with the current version:
library(riskRegression)
library(future)
library(foreach)
library(doFuture)
library(survival)
# Register all available cores for parallel operation
plan(multiprocess, workers = availableCores())
registerDoFuture()
trainSurv <- sampleData(400,outcome="survival")
cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv,
             y=TRUE, x = TRUE)
# Bootstrapping on multiple cores
x1 = Score(list("Cox(X1+X2+X7+X9)"=cox1),
           formula=Surv(time,event)~1,data=trainSurv, times=c(5,8), 
           parallel = "as.registered", split.method="bootcv",B=100)

## End(Not run)
