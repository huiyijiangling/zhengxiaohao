library(Publish)
data(Diabetes)
Diabetes$hyper1 <- factor(1*(Diabetes$bp.1s>140))
## collect odds ratios from three univariate logistic regression analyses 单因素
uni.odds <- glmSeries(hyper1~1,vars=c("chol","hdl","location"),data=Diabetes,family=binomial)
uni.odds
publish(uni.odds)
## control the logistic regression analyses for age and gender 
## but collect only information on the variables in `vars'.
controlled.odds <- glmSeries(hyper1~age+gender,
                             vars=c("chol","hdl","location"),
                             data=Diabetes, family=binomial)
controlled.odds

#这里使用了age+gender作为校正
publish(f)

############################2
## linear regression
data(Diabetes)
f <- glm(bp.1s~AgeGroups+chol+gender+location,data=Diabetes)
rtf <- regressionTable(f,factor.reference = "inline")
plot(rtf,cex=1.3)
## logistic regression
data(Diabetes)
f <- glm(I(BMI>25)~bp.1s+AgeGroups+chol+gender+location,data=Diabetes,family=binomial)
rtf <- regressionTable(f,factor.reference = "inline")
plot(rtf,cex=1.3)
rtf <- regressionTable(f,factor.reference = "extraline")
plot(rtf,cex=1.3)

## Poisson regression
data(trace)
fit <- glm(dead ~ smoking+ sex+ age+Time+offset(log(ObsTime)), family = poisson,data=trace)
rtab <- regressionTable(fit,factor.reference = "inline")
plot(rtab,xlim=c(0.85,1.15),cex=1.8,xaxis.cex=1.5)
## Cox regression
library(survival)
data(pbc)
coxfit <- coxph(Surv(time,status!=0)~age+log(bili)+log(albumin)+factor(edema)+sex,data=pbc)
pubcox <- publish(coxfit)
plot(pubcox,cex=1.5,xratio=c(0.4,0.2))

##########################3
data(Diabetes)
## Linear regression
f = glm(bp.2s~frame+gender+age,data=Diabetes)
publish(f)
publish(f,factor.reference="inline")
publish(f,pvalue.stars=TRUE)
publish(f,ci.format="(l,u)")
### interaction
fit = glm(bp.2s~frame+gender*age,data=Diabetes)
summary(fit)
publish(fit)
Fit = glm(bp.2s~frame*gender+age,data=Diabetes)
publish(Fit)

## Logistic regression
Diabetes$hyper1 <- factor(1*(Diabetes$bp.1s>140))
lrfit <- glm(hyper1~frame+gender+age,data=Diabetes,family=binomial)
publish(lrfit)
### interaction
lrfit1 <- glm(hyper1~frame+gender*age,data=Diabetes,family=binomial)
publish(lrfit1)
lrfit2 <- glm(hyper1~frame*gender+age,data=Diabetes,family=binomial)
publish(lrfit2)
## Poisson regression
data(trace)
trace <- Units(trace,list("age"="years"))
fit <- glm(dead ~ smoking+sex+age+Time+offset(log(ObsTime)), family="poisson",data=trace)
rtf <- regressionTable(fit,factor.reference = "inline")
summary(rtf)
publish(fit)

##############################4
data(Diabetes)
publish(t.test(bp.2s~gender,data=Diabetes))
publish(wilcox.test(bp.2s~gender,data=Diabetes))
publish(with(Diabetes,t.test(bp.2s,bp.1s,paired=TRUE)))
publish(with(Diabetes,wilcox.test(bp.2s,bp.1s,paired=TRUE)))
#############################5
## Not run: 
if (requireNamespace("riskRegression",quietly=TRUE)
    & requireNamespace("mitools",quietly=TRUE)
    & requireNamespace("smcfcs",quietly=TRUE)){
  library(riskRegression)
  library(mitools)
  library(smcfcs)
  ## continuous outcome: linear regression
  # lava some data with missing values
  set.seed(7)
  d=sampleData(78)
  ## generate missing values
  d[X1==1,X6:=NA] 
  d[X2==1,X3:=NA]
  d=d[,.(X8,X4,X3,X6,X7)]
  sapply(d,function(x)sum(is.na(x)))
  
  # multiple imputation (should set m to a large value)
  
  set.seed(17)
  f= smcfcs(d,smtype="lm",
            smformula=X8~X4+X3+X6+X7,
            method=c("","","logreg","norm",""),m=3)
  ccfit=lm(X8~X4+X3+X6+X7,data=d)
  mifit=MIcombine(with(imputationList(f$impDatasets),
                       lm(X8~X4+X3+X6+X7)))
  publish(mifit,fit=ccfit,data=d)
  publish(ccfit)
  
  ## binary outcome
  # lava some data with missing values
  set.seed(7)
  db=sampleData(78,outcome="binary")
  ## generate missing values
  db[X1==1,X6:=NA] 
  db[X2==1,X3:=NA]
  db=db[,.(Y,X4,X3,X6,X7)]
  sapply(db,function(x)sum(is.na(x)))
  
  # multiple imputation (should set m to a large value)
  set.seed(17)
  fb= smcfcs(db,smtype="logistic",
             smformula=Y~X4+X3+X6+X7,
             method=c("","","logreg","norm",""),m=2)
  ccfit=glm(Y~X4+X3+X6+X7,family="binomial",data=db)
  mifit=MIcombine(with(imputationList(fb$impDatasets),
                       glm(Y~X4+X3+X6+X7,family="binomial")))
  publish(mifit,fit=ccfit)
  publish(ccfit)
  
  ## survival: Cox regression
  library(survival)
  # lava some data with missing values
  set.seed(7)
  ds=sampleData(78,outcome="survival")
  ## generate missing values
  ds[X5==1,X6:=NA] 
  ds[X2==1,X3:=NA]
  ds=ds[,.(time,event,X4,X3,X6,X7)]
  sapply(ds,function(x)sum(is.na(x)))
  
  set.seed(17)
  fs= smcfcs(ds,smtype="coxph",
             smformula="Surv(time,event)~X4+X3+X6+X7",
             method=c("","","","logreg","norm",""),m=2)
  ccfit=coxph(Surv(time,event)~X4+X3+X6+X7,data=ds)
  mifit=MIcombine(with(imputationList(fs$impDatasets),
                       coxph(Surv(time,event)~X4+X3+X6+X7)))
  publish(mifit,fit=ccfit,data=ds)
  publish(ccfit)
  
  ## competing risks: Cause-specific Cox regression 
  library(survival)
  # lava some data with missing values
  set.seed(7)
  dcr=sampleData(78,outcome="competing.risks")
  ## generate missing values
  dcr[X5==1,X6:=NA] 
  dcr[X2==1,X3:=NA]
  dcr=dcr[,.(time,event,X4,X3,X6,X7)]
  sapply(dcr,function(x)sum(is.na(x)))
  
  set.seed(17)
  fcr= smcfcs(dcr,smtype="compet",
              smformula=c("Surv(time,event==1)~X4+X3+X6+X7",
                          "Surv(time,event==2)~X4+X3+X6+X7"),
              method=c("","","","logreg","norm",""),m=2)
  ## cause 2 
  ccfit2=coxph(Surv(time,event==2)~X4+X3+X6+X7,data=dcr)
  mifit2=MIcombine(with(imputationList(fcr$impDatasets),
                        coxph(Surv(time,event==2)~X4+X3+X6+X7)))
  publish(mifit2,fit=ccfit2,data=dcr)
  publish(ccfit2)
}

## End(Not run) 
#############################7
data(Diabetes)
f <- glm(bp.1s~age+chol+gender+location,data=Diabetes)
publish(summary(aov(f)),digits=c(1,2))
###############################8

# linear regression
data(Diabetes)
f1 <- glm(bp.1s~age+gender+frame+chol,data=Diabetes)
summary(regressionTable(f1))
summary(regressionTable(f1,units=list("chol"="mmol/L","age"="years")))
## with interaction
f2 <- glm(bp.1s~age*gender+frame+chol,data=Diabetes)
summary(regressionTable(f2))
#Add reference values
summary(regressionTable(f2))
f3 <- glm(bp.1s~age+gender*frame+chol,data=Diabetes)
publish(f3)
regressionTable(f3)
# logistic regression
Diabetes$hyp1 <- factor(1*(Diabetes$bp.1s>140))
l1 <- glm(hyp1~age+gender+frame+chol,data=Diabetes,family="binomial")
regressionTable(l1)
publish(l1)
plot(regressionTable(l1))
## with interaction
l2 <- glm(hyp1~age+gender+frame*chol,data=Diabetes,family="binomial")
regressionTable(l2)
l3 <- glm(hyp1~age*gender+frame*chol,data=Diabetes,family="binomial")
regressionTable(l3)
# Cox regression
library(survival)
data(pbc)
pbc$edema <- factor(pbc$edema,levels=c("0","0.5","1"),labels=c("0","0.5","1"))
c1 <- coxph(Surv(time,status!=0)~log(bili)+age+protime+sex+edema,data=pbc)
regressionTable(c1)
# with interaction
c2 <- coxph(Surv(time,status!=0)~log(bili)+age+protime*sex+edema,data=pbc)
regressionTable(c2)
c3 <- coxph(Surv(time,status!=0)~edema*log(bili)+age+protime+sex+edema+edema:sex,data=pbc)
regressionTable(c3)
if (requireNamespace("nlme",quietly=TRUE)){
  ## gls regression
  library(lava)
  library(nlme)
  m <- lvm(Y ~ X1 + gender + group + Interaction)
  distribution(m, ~gender) <- binomial.lvm()
  distribution(m, ~group) <- binomial.lvm(size = 2)
  constrain(m, Interaction ~ gender + group) <- function(x){x[,1]*x[,2]}
  d <- sim(m, 1e2)
  d$gender <- factor(d$gender, labels = letters[1:2])
  d$group <- factor(d$group)
  e.gls <- gls(Y ~ X1 + gender*group, data = d,
               weights = varIdent(form = ~1|group))
  regressionTable(e.gls)
  summary(regressionTable(e.gls))
}
#######################################################9
#load libraries
library(data.table)
library(Publish)
library(survival)
data(traceR) #get dataframe traceR
data.table::setDT(traceR)
traceR[,':='(wmi2=factor(wallMotionIndex<0.9,levels=c(TRUE,FALSE),
                         labels=c("bad","good")),
             abd2=factor(abdominalCircumference<95, levels=c(TRUE,FALSE),
                         labels=c("slim","fat")))]
traceR[,sex:=as.factor(sex)] # all subgroup variables needs to be factor
traceR[observationTime==0,observationTime:=1]
# remove missing covariate values
traceR=na.omit(traceR)
# univariate analysis of smoking in subgroups of age and sex
# Main regression analysis is a simple/univariate Cox regression
fit_cox <- coxph(Surv(observationTime,dead)~treatment,data=traceR)
sub_cox <- subgroupAnalysis(fit_cox,traceR,treatment="treatment",
                            subgroups=c("smoking","sex","wmi2","abd2"))
sub_cox
# to see how the results are obtained consider the variable: smoking
fit_cox_smoke <- coxph(Surv(observationTime,dead)~treatment*smoking,data=traceR)
# the last three rows of the following output:
publish(fit_cox_smoke)
# are included in the first 3 rows of the result of the sub group analysis:
sub_cox[1:3,]
# the p-value is obtained as:
fit_cox_smoke_add <- coxph(Surv(observationTime,dead)~treatment+smoking,data=traceR)
anova(fit_cox_smoke_add,fit_cox_smoke,test="Chisq")
# Note that a real subgroup analysis would be to subset the data
fit_cox1a <- coxph(Surv(observationTime,dead)~treatment,data=traceR[smoking=="never"])
fit_cox1b <- coxph(Surv(observationTime,dead)~treatment,data=traceR[smoking=="current"])
fit_cox1c <- coxph(Surv(observationTime,dead)~treatment,data=traceR[smoking=="prior"])
## when the main analysis is already adjusted
fit_cox_adj <- coxph(Surv(observationTime,dead)~treatment+smoking+sex+wmi2+abd2,
                     data=traceR)
sub_cox_adj <- subgroupAnalysis(fit_cox_adj,traceR,treatment="treatment",
                                subgroups=c("smoking","sex","wmi2","abd2")) # subgroups as character string
sub_cox_adj
# When both start and end are in the Surv statement:
traceR[,null:=0]
fit_cox2 <- coxph(Surv(null,observationTime,dead)~treatment+smoking+sex+wmi2+abd2,data=traceR)
summary(regressionTable(fit_cox))
sub_cox2 <- subgroupAnalysis(fit_cox2,traceR,treatment="treatment",
                             subgroups=c("smoking","sex","wmi2","abd2"))
# Analysis with Poisson - and the unrealistic assumption of constant hazard
# and adjusted for age in all subgroups
fit_p <- glm(dead~treatment+age+offset(log(observationTime)),family="poisson",
             data=traceR)
sub_pois <- subgroupAnalysis(fit_p,traceR,treatment="treatment",
                             subgroups=~smoking+sex+wmi2+abd2)
# Analysis with logistic regression - and very wrongly ignoring censoring
fit_log <- glm(dead~treatment+age,family="binomial",data=traceR)
sub_log <- subgroupAnalysis(fit_log,traceR,treatment="treatment",
                            subgroups=~smoking+sex+wmi2+abd2, factor.reference="inline")
###################################################10
library(survival)
data(pbc)
pbc$edema <- factor(pbc$edema,levels=c("0","0.5","1"),labels=c("0","0.5","1"))
fit = coxph(Surv(time,status!=0)~age+sex+edema+log(bili)+log(albumin)+log(protime),
            data=pbc)
u=summary(regressionTable(fit))
u$regressionTable
u$rawTable
summary(regressionTable(fit),handler="prettyNum")
summary(regressionTable(fit),handler="format")
summary(regressionTable(fit),handler="sprintf",digits=c(2,2),pValue.stars=TRUE)
summary(regressionTable(fit),handler="sprintf",digits=c(2,2),pValue.stars=TRUE,ci.format="(l,u)")
################################################11

data(Diabetes)
library(data.table)
univariateTable(~age,data=Diabetes)
univariateTable(~gender,data=Diabetes)
univariateTable(~age+gender+ height+weight,data=Diabetes)
## same thing but less typing
utable(~age+gender+ height+weight,data=Diabetes)
## summary by location:
univariateTable(location~Q(age)+gender+height+weight,data=Diabetes)
## continuous variables marked with Q() are (by default) summarized
## with median (IQR) and kruskal.test (with two groups equivalent to wilcox.test)
## variables not marked with Q() are (by default) summarized
## with mean (sd) and anova.glm(...,test="Chisq")
## the p-value of anova(glm()) with only two groups is similar
## but not exactly equal to that of a t.test
## categorical variables are (by default) summarized by count
## (percent) and chi-square tests (\code{chisq.test}). When \code{compare.groups ='logistic'}
## anova(glm(...,family=binomial,test="Chisq")) is used to calculate p-values.
## export result to csv
table1 = summary(univariateTable(location~age+gender+height+weight,data=Diabetes),
                 show.pvalues=FALSE)
# write.csv(table1,file="~/table1.csv",rownames=FALSE)
## change labels and values
utable(location~age+gender+height+weight,data=Diabetes,
       age="Age (years)",gender="Sex",
       gender.female="Female",
       gender.male="Male",
       height="Body height (inches)",
       weight="Body weight (pounds)")
## Use quantiles and rank tests for some variables and mean and standard deviation for others
univariateTable(gender~Q(age)+location+Q(BMI)+height+weight,
                data=Diabetes)
## Factor with more than 2 levels
Diabetes$AgeGroups <- cut(Diabetes$age,
                          c(19,29,39,49,59,69,92),
                          include.lowest=TRUE)
univariateTable(location~AgeGroups+gender+height+weight,
                data=Diabetes)
summary(univariateTable(location~AgeGroups+gender+height+weight,
                data=Diabetes))
summary(univariateTable(location~AgeGroups+gender+height+weight,
                        data=Diabetes),show.missing = "never")
## Row percent
univariateTable(location~gender+age+AgeGroups,
                data=Diabetes,
                column.percent=FALSE)
## change of frequency format
univariateTable(location~gender+age+AgeGroups,
                data=Diabetes,
                column.percent=FALSE,
                freq.format="percent(x) (n=count(x))")
## changing Labels
u <- univariateTable(location~gender+AgeGroups+ height + weight,
                     data=Diabetes,
                     column.percent=TRUE,
                     freq.format="count(x) (percent(x))")
summary(u,"AgeGroups"="Age (years)","height"="Height (inches)")
## more than two groups
Diabetes$frame=factor(Diabetes$frame,levels=c("small","medium","large"))
univariateTable(frame~gender+BMI+age,data=Diabetes)
Diabetes$sex=as.numeric(Diabetes$gender)
univariateTable(frame~sex+gender+BMI+age,
                data=Diabetes,freq.format="count(x) (percent(x))")
## multiple summary formats
## suppose we want for some reason mean (range) for age
## and median (range) for BMI.
## method 1:
univariateTable(frame~Q(age)+BMI,
                data=Diabetes,
                Q.format="mean(x) (range(x))",
                summary.format="median(x) (range(x))")
## method 2:
u1 <- summary(univariateTable(frame~age,
                              data=na.omit(Diabetes),
                              summary.format="mean(x) (range(x))"))
u2 <- summary(univariateTable(frame~BMI,
                              data=na.omit(Diabetes),
                              summary.format="median(x) (range(x))"))
publish(rbind(u1,u2),digits=2)
## Large number format (big.mark)
Diabetes$AGE <- 1000*Diabetes$age
u3 <- summary(univariateTable(frame~AGE,
                              data=Diabetes,big.mark="'"))
########################
