library(survival)
library(survminer)
NCOA1_G_S <- read.csv("Survial.CSV")
group <- ifelse(NCOA1_G_S $NCOA1>median(NCOA1_G_S$NCOA1)&NCOA1_G_S $GLI2>median(NCOA1_G_S$GLI2),'Nh-Gh',
                ifelse(NCOA1_G_S $NCOA1<median(NCOA1_G_S$NCOA1)&NCOA1_G_S$GLI2<median(NCOA1_G_S$GLI2),'NL-GL',
                       ifelse(NCOA1_G_S$NCOA1>median(NCOA1_G_S$NCOA1)&NCOA1_G_S$GLI2<median(NCOA1_G_S$GLI2),'Nh-Gl','NL-Gh')))

sfit <- survfit(Surv(Time, Event)~group, data=NCOA1_G_S)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)



MTC是一个非常具有临床意义的选题，与其他预后较好的甲状腺癌相比，尽管必须