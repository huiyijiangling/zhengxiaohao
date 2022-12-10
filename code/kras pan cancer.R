library(dplyr)
library(survival)
library(survminer)
library(ggtext)

krasall_muttable = read.table( r'(C:\Users\zxh\Desktop\R\paad-tcga-gtex\kras mut cohort\allkrasmut\alterations_across_samples (2).tsv)',sep = '\t',header = T)
krasall_surtable = read.table( r'(C:\Users\zxh\Desktop\R\paad-tcga-gtex\kras mut cohort\allkrasmut\Overall_ .txt)',sep = '\t',header = T,skip = 3)
krasall_muttable%>%group_by(KRAS..MUT) %>% mutate(count = n()) ->krasall_muttable

krasall_muttable$DriverPassenger=ifelse(grepl("driver",krasall_muttable$KRAS..MUT,ignore.case = T),"Driver mutation","Passenger mutation")
krasall_muttable$complex_mut=ifelse(grepl("\\,|\\*|del|ins|splice|dup",krasall_muttable$KRAS..MUT,ignore.case = T),"Complex","Substitution")
dd=subset(krasall_muttable,complex_mut=="Substitution")
table(dd$KRAS..MUT)
if(F){
krasall_muttable=subset(krasall_muttable,count>400)#snp具体
krasall_muttable=subset(krasall_muttable,krasall_muttable$KRAS..AMP !="not profiled")#cnv
}

krasall_muttable[krasall_muttable$KRAS..AMP =="AMP (driver), AMP","KRAS..AMP"]="AMP (driver)"

krasall_surtable=subset(krasall_surtable,!is.na(krasall_surtable$Status))

res.cat.va=merge(krasall_muttable,krasall_surtable,by.x="Patient.ID",by.y = "Case.ID",all.x = F,all.y = F)

table(res.cat.va$Status,useNA = "ifany")
table(krasall_muttable$KRAS..MUT)
table(krasall_muttable$count)
table(krasall_muttable$KRAS..AMP)

res.cat.va$event=as.numeric(ifelse(res.cat.va$Status=="censored",0,1))
res.cat.va$time=as.numeric(res.cat.va$Time..months.)

fit <- survival::survfit(survival::Surv(time, event) ~`complex_mut`, data = res.cat.va)
x="complex_mut"
#DriverPassenger
#`KRAS..MUT`
#`KRAS..AMP`
# s=as.formula(paste0("survival::Surv(time, event) ~",as.character(x)))
# model <- survival::coxph(formula=s, data =dat)
# # summary(model,data=dat)
# # options(scipen=1)
# new_dat=dat
# new_dat$event=as.numeric(new_dat$event)
# new_dat$time=as.numeric(new_dat$time)
# fp <- predict(model,new_dat)
# library(cowplot)
# library(pheatmap)
# risk=new_dat
# risk$fp=fp
# res.cut.va <- surv_cutpoint(risk, time = "time", event = "event",minprop=0.2,
#                             variables =c("fp"))
# summary(res.cut.va)
# res.cat.va <- survminer::surv_categorize(res.cut.va)
# fit <- survival::survfit(survival::Surv(time, event) ~fp, data = res.cat.va)

# pdf(file=paste0("./piliangshengcuntu/","Figure.",x,"_best_survplot.pdf"),height = 6,width = 6,onefile = T)
edr1=survminer::ggsurvplot(fit,
                           data = res.cat.va,
                           legend.title = x,#as.character(choose_gene[choose_gene$Gene%in% as.character(x),"GeneName"]),
                           # legend.labs = c( "High","Low"),
                           risk.table = T,#循环时不能用否则会报错，单个人可以T
                           conf.int = TRUE,
                           pval=TRUE,
                           tables.height = 0.2,
                           tables.theme = theme_cleantable(),
                           # palette = c("#E7B800","#2E9FDF"),
                           ggtheme = theme_bw()
)
ggsurv=edr1
g2 <- ggplotGrob(ggsurv$plot )
g3 <- ggplotGrob(ggsurv$table)
min_ncol <- min(ncol(g2), ncol(g3))
g <- gridExtra::gtable_rbind(g2[, 1:min_ncol], g3[, 1:min_ncol], size="last")
g$widths <- grid::unit.pmax(g2$widths, g3$widths)
grid::grid.newpage()
grid::grid.draw(g)
# ggsave(filename = paste0("Figure.",x,"_best_survplot.pdf"),
#        plot=print(edr1, newpage = FALSE),
#        path = paste0("./piliangshengcuntu/",stringr::str_split(multiplesurvival,"\\/",simplify=T)[,8],"/bestcut"),
#        device = 'pdf',
#        width = 8,height=8,dpi = 600)
ggsave(filename = paste0("Figure. KRAS",x," amp or not.pdf"),
       g,
       path = "./",
       device = 'pdf',
       width = 8,height=8,dpi = 600)
dev.off()
# ggsave(filename = paste0("Figure.",x,"_best_survplot.pdf"),
#        g,
#        path = "./",
#        device = 'pdf',
#        width = 8,height=8,dpi = 600)



edr1=survminer::ggsurvplot(fit,
                           data = res.cat.va,
                           legend.title = x,#as.character(choose_gene[choose_gene$Gene%in% as.character(x),"GeneName"]),
                           # legend.labs = c( "High","Low"),
                           risk.table = F,#循环时不能用否则会报错，单个人可以T
                           conf.int = TRUE,
                           pval=TRUE,
                           tables.height = 0.2,
                           tables.theme = theme_cleantable(),
                           # palette = c("#E7B800","#2E9FDF"),
                           ggtheme = theme_bw()
)
edr1
ggsave(filename = paste0("Figure.",x,"fff_best_survplot.pdf"),
       path = "./",
       device = 'pdf',
       width = 8,height=8,dpi = 600)
dev.off()
