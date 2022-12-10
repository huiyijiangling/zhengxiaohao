library(survival)
library(survminer)

res.cut <- surv_cutpoint(dat, time = "time", event = "event",
                         variables =c('TTF2','TXNDC9','ZNF562','VEZT','YY1','MED4'))
summary(res.cut)
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
ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE)
