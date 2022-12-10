###############################################aftersurgery
JDW_replustime_huayan_aftersurgery=subset(JDW_replustime_huayan,as.Date(JDW_replustime_huayan$标本采集时间,tz=Sys.timezone(location = TRUE))>as.Date(JDW_replustime_huayan$ssrq,tz=Sys.timezone(location = TRUE)))#忽略了多个当天的抽的
# JDW_replustime_huayan_aftersurgery=subset(JDW_replustime_huayan,as.Date(JDW_replustime_huayan$标本采集时间,tz=Sys.timezone(location = TRUE))<(as_datetime(JDW_replustime_huayan$ssrq,tz=Sys.timezone(location = TRUE))+8.5*60*60))#忽略了多个当天的抽的
length(unique(JDW_replustime_huayan_aftersurgery$病案号))
JDW_replustime_huayan_dangtianchouxue=subset(JDW_replustime_huayan,as.Date(JDW_replustime_huayan$标本采集时间,tz=Sys.timezone(location = TRUE))==as.Date(JDW_replustime_huayan$ssrq,tz=Sys.timezone(location = TRUE)))
# 519 人有化验
JDW_replustime_huayan_aftersurgery_nearest=arrange(JDW_replustime_huayan_aftersurgery,desc(JDW_replustime_huayan_aftersurgery$标本采集时间))
JDW_replustime_huayan_aftersurgery_nearest=JDW_replustime_huayan_aftersurgery_nearest[!duplicated(JDW_replustime_huayan_aftersurgery_nearest[,c("项目名称","病案号")]),]
###############


box=JDW_replustime_huayan_aftersurgery_nearest
# box$项目结果=log2(box$项目结果+1)
box=pivot_wider(
  box,
  id_cols = c("病案号","zyts","sqts", "shts" ,"qzryts", "qzssts","yqts", "dsts"),#,"neoadjuvant","After2020"
  names_from = "项目名称",
  names_prefix = "",
  # names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "check_unique",
  values_from = "项目结果",
  values_fill = NULL,
  values_fn = NULL)
box=as.data.frame(box)
rownames(box) = box$病案号
box$病案号=NULL
library(corrplot)
box[]=lapply(box, as.numeric)
res <- list()
res[[1]]=box
res[[2]]=box[rownames(box)%in%JDW_replustime[JDW_replustime$neoadjuvant==0,]$bah,]
res[[3]]=box[rownames(box)%in%JDW_replustime[JDW_replustime$neoadjuvant==1,]$bah,]
res[[4]]=box[rownames(box)%in%JDW_replustime[JDW_replustime$After2020==0,]$bah,]
res[[5]]=box[rownames(box)%in%JDW_replustime[JDW_replustime$After2020==1,]$bah,]
pp11=foreach(i=1:5)%do%{
  pdf(paste0("CoR",i,".pdf"),width = 10,height=10)
  # box=JDW_replustime_huayan_aftersurgery_nearest
  # box$项目结果=log2(box$项目结果+1)
  res1 <- cor.mtest(res[[i]], conf.level = .95, method = "spearman")
  pp1<- cor(res[[i]], method ="spearman",use="complete.obs") %>%  corrplot(type="lower",method = "square",
                                                                           insig="label_sig",p.mat = res1$p,outline = "white",
                                                                           sig.level=c(.05),  ,#.001,.01,
                                                                           pch.cex=.9,pch.col="green",tl.srt=45,tl.col="black")
  dev.off()
  pp11=pp1$corrPos
  pp11=pp11[which(pp11$xName%in%c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "yqts", "dsts")|pp11$yName%in%c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "yqts", "dsts")),]
  pp11=subset(pp11,xName!=yName)
  pp11=subset(pp11,!pp11$yName%in%c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "yqts", "dsts"))
  pp11=unite(pp11,col = "ts_husyan",xName,yName,remove = F)
  openxlsx::write.xlsx(pp11,paste0("天数与术前化验相关性",i,".xlsx"))
  pp11_cor=pivot_wider(
    pp11,
    id_cols = c("yName"),
    names_from = "xName",
    names_prefix = "",
    # names_sep = "_",
    names_glue = NULL,
    names_sort = FALSE,
    names_repair = "check_unique",
    values_from = "corr",
    values_fill = NULL,
    values_fn = NULL)
  pp11_cor=data.frame(pp11_cor,row.names = "yName")
  pp11_pvalue=pivot_wider(
    pp11,
    id_cols = c("yName"),
    names_from = "xName",
    names_prefix = "",
    # names_sep = "_",
    names_glue = NULL,
    names_sort = FALSE,
    names_repair = "check_unique",
    values_from = "p.value",
    values_fill = NULL,
    values_fn = NULL)
  pp11_pvalue=data.frame(pp11_pvalue,row.names = "yName")
  pp11=pp11[which(pp11$p.value<0.05),]
  return(list(pp11,pp11_cor,pp11_pvalue))
}

ts_husyan_name <- list()
ts_husyan_sign <- list()
common_cor <- list()
for (i in 1:5) {
  common_cor[[i]]=pp11[[i]][[1]]
  ts_husyan_name[[i]]=common_cor[[i]]$ts_husyan  
}
output_huayanzhibiao=unique(rbindlist(common_cor)$yName)
output_huayanzhibiao=output_huayanzhibiao[output_huayanzhibiao%in%unique(c(common_cor[[2]]$yName,common_cor[[3]]$yName))]
dput(output_huayanzhibiao)
output_huayanzhibiao=output_huayanzhibiao[output_huayanzhibiao%in%c("ALT", "AST", "CK-MB", "CRE", "G", "GGT", "GLU","HBDH", "PALB", "TP", "EOS#", "EOS%", "MONO#", "MONO%","MPV", "P-LCR", "PDW", "SOD", "IBIL","TRANSFE", "HCT",  "ALB", "ALP", "IGG", "LDH","TBA", "MCH", "MCV", "RBC", "RDW-CV", "RDW-SD","WBC", "HB", "UREA", "DBIL", "CK")]
# "PT(A)", "PT(R)", "PT(S)","APTT", "AT-III", "D-D", "FDP", "PT(INR)", "TT",
# "ADA",
writexl::write_xlsx(common_cor,path = "common_cor.xlsx")
########    
ts_husyan_name12=Reduce(intersect,ts_husyan_name[c(1,2)])
ts_husyan_sign <- list()
for (i in c(1,2)) {
  ts_husyan_sign[[i]]=common_cor[[i]][common_cor[[i]]$ts_husyan%in%ts_husyan_name12,]
  ts_husyan_sign[[i]]=ts_husyan_sign[[i]][match(ts_husyan_sign[[i]]$ts_husyan,ts_husyan_name12),]
}
# ts_husyan_sign[[2]]=NULL
ts_husyan_sign=as.data.frame(ts_husyan_sign)
ts_husyan_sign22=subset(ts_husyan_sign,!ts_husyan_sign$yName%in%c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "yqts", "dsts"))
ts_husyan_sign22=subset(ts_husyan_sign22,abs(sign(ts_husyan_sign22$corr)+sign(ts_husyan_sign22$corr.1))==2)

# ts_husyan_sign22=subset(ts_husyan_sign22,abs(ts_husyan_sign22$corr)>0.2&abs(ts_husyan_sign22$corr.1)>0.2)
openxlsx::write.xlsx(ts_husyan_sign22,file ="ts_husyan_sign22.xlsx")
###########################
ts_husyan_name13=Reduce(intersect,ts_husyan_name[c(1,3)])
ts_husyan_sign <- list()
for (i in c(1,3)) {
  ts_husyan_sign[[i]]=common_cor[[i]][common_cor[[i]]$ts_husyan%in%ts_husyan_name13,]
  ts_husyan_sign[[i]]=ts_husyan_sign[[i]][match(ts_husyan_sign[[i]]$ts_husyan,ts_husyan_name13),]
}
ts_husyan_sign[[2]]=NULL
ts_husyan_sign=as.data.frame(ts_husyan_sign)
ts_husyan_sign33=subset(ts_husyan_sign,!ts_husyan_sign$yName%in%c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "yqts", "dsts"))
ts_husyan_sign33=subset(ts_husyan_sign33,abs(sign(ts_husyan_sign33$corr)+sign(ts_husyan_sign33$corr.1))==2)
# ts_husyan_sign33=subset(ts_husyan_sign33,abs(ts_husyan_sign33$corr)>0.2&abs(ts_husyan_sign33$corr.1)>0.2)
openxlsx::write.xlsx(ts_husyan_sign33,file ="ts_husyan_sign33.xlsx")
######################

# ts_husyan_sign22=subset(ts_husyan_sign22,abs(sign(ts_husyan_sign22$corr)+sign(ts_husyan_sign22$corr.1)+sign(ts_husyan_sign22$corr.2))==3)
if(T){
  box=JDW_replustime_huayan_aftersurgery_nearest
  # box$项目结果=log2(box$项目结果+1)
  box=tidyr::pivot_longer(
    data=box,
    cols=c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "yqts", "dsts"),#everthing
    names_to = c("DAYtype"),
    names_prefix = NULL,
    # names_sep = "_",
    names_pattern = NULL,
    # names_ptypes = list(),
    # names_transform = list(),
    names_repair = "check_unique",
    values_to = "Daytime",
    values_drop_na = T,
    # values_ptypes = list(),
    # values_transform = list()
  )
  colnames(box)[duplicated(colnames(box))]
  pdf(file =paste0("Figure A. 化验全因相关.pdf"),height = 200,width = 30)
  box$Daytime=as.numeric(box$Daytime)
  # box$Daytime=log2(box$Daytime+1)
  p1=ggplot(box, aes(x = 项目结果, y = Daytime)) + 
    ylab("")+xlab("")+
    geom_point(shape = 21, colour = "#4682B4", fill = "#87CEFA", size = 3, stroke = .5,alpha=0.8)+ geom_smooth(method="glm",formula = y ~ x,linetype=2,color="#6495ED",fill="#D3D3D3") + theme_bw()+stat_cor(method = 'spearman', aes(),na.rm = T,size = 10)+facet_grid(rows=vars(box$项目名称),cols=vars(box$DAYtype),margins = F)#,space="free_y"
  # p2=ggExtra::ggMarginal(p1, type = "density", xparams = list(fill = "#FFE4B5"),yparams = list(fill = "#90EE90"))
  p1
  dev.off()
}
JDW_replustime_huayan_aftersurgery_nearest <- JDW_replustime_huayan_aftersurgery_nearest %>%
  mutate(across(where(is.difftime),as.numeric)) 

quantile(JDW_replustime_huayan_aftersurgery_nearest$sqts)
JDW_replustime_huayan_aftersurgery_nearest$sqts1 = cut(JDW_replustime_huayan_aftersurgery_nearest$sqts,breaks = c(-Inf,4,Inf),right =T,ordered_result = T,labels = c("≤4",">4"))
quantile(JDW_replustime_huayan_aftersurgery_nearest$shts)
JDW_replustime_huayan_aftersurgery_nearest$shts1 = cut(JDW_replustime_huayan_aftersurgery_nearest$shts,breaks = c(-Inf,10,Inf),right =T,ordered_result = T,labels = c("≤10",">10"))
quantile(JDW_replustime_huayan_aftersurgery_nearest$zyts)
JDW_replustime_huayan_aftersurgery_nearest$zyts1 = cut(JDW_replustime_huayan_aftersurgery_nearest$zyts,breaks = c(-Inf,14,Inf),right = T,ordered_result = T,labels = c("≤14",">14"))
quantile(JDW_replustime_huayan_aftersurgery_nearest$yqts)
JDW_replustime_huayan_aftersurgery_nearest$yqts1 = cut(JDW_replustime_huayan_aftersurgery_nearest$yqts,breaks = c(-Inf,14,Inf),right =T,ordered_result = T,labels = c("≤14",">14"))
JDW_replustime_huayan_aftersurgery_nearest$yqts2 = cut(JDW_replustime_huayan_aftersurgery_nearest$yqts,breaks = c(-Inf,30,Inf),right =T,ordered_result = T,labels = c("≤30",">30"))                                                                                             
quantile(JDW_replustime_huayan_aftersurgery_nearest$qzryts)
JDW_replustime_huayan_aftersurgery_nearest$qzryts1 = cut(JDW_replustime_huayan_aftersurgery_nearest$qzryts,breaks = c(-Inf,14,Inf),right =T,ordered_result = T,labels = c("≤14",">14"))
quantile(JDW_replustime_huayan_aftersurgery_nearest$qzssts)
JDW_replustime_huayan_aftersurgery_nearest$qzssts1 = cut(JDW_replustime_huayan_aftersurgery_nearest$qzssts,breaks = c(-Inf,21,Inf),right =T,ordered_result = T,labels = c("≤21",">21"))
quantile(JDW_replustime_huayan_aftersurgery_nearest$dsts,probs = seq(0, 1, 0.1))
quantile(JDW_replustime_huayan_aftersurgery_nearest[which(JDW_replustime_huayan_aftersurgery_nearest$neoadjuvant==1),]$dsts,probs = seq(0, 1, 0.1))
quantile(JDW_replustime_huayan_aftersurgery_nearest[which(JDW_replustime_huayan_aftersurgery_nearest$neoadjuvant==0),]$dsts,probs = seq(0, 1, 0.1))
JDW_replustime_huayan_aftersurgery_nearest$dsts1 = cut(JDW_replustime_huayan_aftersurgery_nearest$dsts,breaks = c(-Inf,14,Inf),right =T,ordered_result = T,labels =c("≤14",">14"))
JDW_replustime_huayan_aftersurgery_nearest$dsts2 = cut(JDW_replustime_huayan_aftersurgery_nearest$dsts,breaks = c(-Inf,30,Inf),right =T,ordered_result = T,labels =c("≤30",">30"))
JDW_replustime_huayan_aftersurgery_nearest$dsts3=NA
JDW_replustime_huayan_aftersurgery_nearest[which(JDW_replustime_huayan_aftersurgery_nearest$neoadjuvant==1),]$dsts3=ifelse(JDW_replustime_huayan_aftersurgery_nearest[which(JDW_replustime_huayan_aftersurgery_nearest$neoadjuvant==1),]$dsts>30,1,ifelse(JDW_replustime_huayan_aftersurgery_nearest[which(JDW_replustime_huayan_aftersurgery_nearest$neoadjuvant==1),]$dsts<=30,0,"wrong"))

JDW_replustime_huayan_aftersurgery_nearest[which(JDW_replustime_huayan_aftersurgery_nearest$neoadjuvant==0),]$dsts3=ifelse(JDW_replustime_huayan_aftersurgery_nearest[which(JDW_replustime_huayan_aftersurgery_nearest$neoadjuvant==0),]$dsts>14,1,ifelse(JDW_replustime_huayan_aftersurgery_nearest[which(JDW_replustime_huayan_aftersurgery_nearest$neoadjuvant==0),]$dsts<=14,0,"wrong"))
JDW_replustime_huayan_aftersurgery_nearest$dsts3=as.numeric(JDW_replustime_huayan_aftersurgery_nearest$dsts3)
# JDW_replustime_huayan_aftersurgery_nearest$dsts3=as.factor( JDW_replustime_huayan_aftersurgery_nearest$dsts3)

box=JDW_replustime_huayan_aftersurgery_nearest
write.csv(JDW_replustime_huayan_aftersurgery_nearest,file="box.csv")
# box$项目结果=log2(box$项目结果+1)
box=pivot_wider(
  box,
  id_cols = c("病案号","neoadjuvant","Sex", "Age", "TS", "size", "Lauren", "Borrmann", "Diff", "signet","pT","pN","pM","Smoking", "Drinking", "HD", "hypertension", "diabetes","Age","Sex","pT","ssrq","zyts","sqts", "shts" ,"qzryts", "qzssts","yqts", "dsts","zyts1","sqts1", "shts1" ,"qzryts1", "qzssts1",  "yqts1", "dsts1","dsts2","dsts3"),#everthing#,"neoadjuvant","After2020"
  names_from = "项目名称",
  names_prefix = "",
  # names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "check_unique",
  values_from = "项目结果",
  values_fill = NULL,
  values_fn = NULL)
unique(box$病案号)

fit_p <- glm(dsts3~GGT+AST+TT+MCV,family=binomial(),data=box)
fit_p <- glm(dsts3~GGT+AST+TT+MCV,family=binomial(),data=box)
summary(fit_p)
fit_p <- MASS::glm.nb(dsts~GGT+AST+TT+MCV+MPV,data=box)
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



p <- ggpubr::ggboxplot(box, x = "项目结果", y = "Daytime",color = "res",#group="项目名称", 
                       font.label = list(size = 1, color = "black"),outlier.shape = NA#,palette=c(1,2,3),order=c("No changes","Shallow Deletion","Deep Deletion")#"#00AFBB","#E7B800",
)+labs(title="",x=" ", y = "log2 ( Laboratory results + 1)")+stat_compare_means(method = "kruskal.test", label.y =9)+stat_compare_means(aes(color="res"),label.y =seq(5, 9, length.out = 15),method = "wilcox.test", comparisons =combn(1:6, 2, FUN = list))+facet_grid(rows = vars(box$项目名称),space="free_y",margins = F)#+facet_wrap(~项目名称,nrow = 6)#+#tag="ENSG00000000003"
p
#主要问题在2020-08-06 2020-08-07，而前面比较相似
dev.off()
}


JDW_replustime_huayan_farest=arrange(JDW_replustime_huayan,JDW_replustime_huayan$标本采集时间)
JDW_replustime_huayan_farest=JDW_replustime_huayan_farest[!duplicated(JDW_replustime_huayan_farest[,c("项目名称","病案号")]),]
length(unique(JDW_replustime_huayan_farest$病案号))

JDW_replustime_huayan_farest=arrange(JDW_replustime_huayan,JDW_replustime_huayan$标本采集时间)
JDW_replustime_huayan_farest=JDW_replustime_huayan_farest[!duplicated(JDW_replustime_huayan_farest[,c("项目名称","病案号")]),]
length(unique(JDW_replustime_huayan_farest$病案号))
#
JDW_replustime_huayan_nearest=arrange(JDW_replustime_huayan,desc(JDW_replustime_huayan$标本采集时间))
JDW_replustime_huayan_nearest=JDW_replustime_huayan_nearest[!duplicated(JDW_replustime_huayan_nearest[,c("项目名称","病案号")]),]
length(unique(JDW_replustime_huayan_nearest$病案号))
save(JDW_replustime_huayan_nearest,JDW_replustime_huayan_farest,JDW_replustime_huayan,file = "JDW_replustime_huayan.Rdata")
# load(file = "JDW_replustime_huayan.Rdata")


#
JDW_replustime_huayan_nearest=arrange(JDW_replustime_huayan,desc(JDW_replustime_huayan$标本采集时间))
JDW_replustime_huayan_nearest=JDW_replustime_huayan_nearest[!duplicated(JDW_replustime_huayan_nearest[,c("项目名称","病案号")]),]
length(unique(JDW_replustime_huayan_nearest$病案号))
save(JDW_replustime_huayan_nearest,JDW_replustime_huayan_farest,JDW_replustime_huayan,file = "JDW_replustime_huayan.Rdata")
# load(file = "JDW_replustime_huayan.Rdata")
library(tidyr)
JDW_replustime_huayan_nearest_expr=pivot_wider(
  JDW_replustime_huayan_nearest,
  id_cols = "病案号",
  names_from = "项目名称",
  names_prefix = "",
  names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "check_unique",
  values_from = "项目结果",
  # values_fill = NULL,
  values_fn = NULL)
writexl::write_xlsx(JDW_replustime_huayan_nearest_expr,"JDW_replustime_huayan_nearest_expr.xlsx")
# unlist(lapply(JDW_replustime_huayan_nearest_expr,function(x) sum(is.na(x))/length(x)>0.25))
wula=JDW_replustime_huayan_nearest_expr[unlist(lapply(JDW_replustime_huayan_nearest_expr,function(x) sum(is.na(x))/length(x)<0.5))]
# "标本采集时间",

# 基函数：x设置目标变量
library(ggplot2)
ggplot(JDW, aes(x = yqts, fill = NAC)) +
  # 直方图函数：position设置堆积模式为重叠
  geom_histogram(position = "identity", alpha = 0.4)
table(JDW_replustime_huayan_nearest$yqts,JDW_replustime_huayan_nearest$NAC)
boxplot(yqts~NAC,JDW)





#术前的情况 是直接按照我们生产时才看  缺不缺 在上面

#而质控时直接使用即可

JDW_replustime_huayan_nearest_expr=pivot_wider(
  box,
  id_cols = "病案号",
  names_from = "项目名称",
  names_prefix = "",
  names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "check_unique",
  values_from = "项目结果",
  # values_fill = NULL,
  values_fn = NULL)






# 
# shapiro.test#3-5000
# lillie.test#>5000
# library(nortest)
# table(box_surv$cnv3)
# 
# tapply(box$MIR31HG..mRNA.Expression..RSEM..Batch.normalized.from.Illumina.HiSeq_RNASeqV2...log2.value...1..,box$cnv3,lillie.test)#不正态
# bartlett.test(ENSG00000171889~cnv,data=box)#方差齐性
# box$标本采集时间
# box$项目结果=log10(box$项目结果+1)

#
# library(doParallel) #
# cl <- makeCluster(30)
# registerDoParallel(cl)
# foreach(x=choose_gene[[2]],.packages = c(library(survival),library(survminer),library(ggtext),library(cowplot))) %dopar% multiple_box(x) # returns list
# 






if(F){
  pdf('Figure A.Boxpgrgrl.pdf',width = 20,height=400)
  box=resultALL_ref_sva
  box$标本采集时间=quarter(box$标本采集时间,with_year = T)
  if(require('ggpubr')){
    library(ggpubr)
    # google search : ggpubr boxplot add p-value
    # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
    p <- ggboxplot(box, x = "标本采集时间", y = "项目结果",#group="项目名称",# color = "项目名称", 
                   font.label = list(size = 1, color = "black"),outlier.shape = NA#,palette=c(1,2,3),order=c("No changes","Shallow Deletion","Deep Deletion")#"#00AFBB","#E7B800",
    )+#add = "jitter"
      # scale_color_manual(values=c(1,2,3,4))  +#"#E69F00", "#56B4E9"
      #    coord_cartesian(ylim = c(0.0000, 30))+
      # theme(legend.position="none")+
      labs(title="",x="MIR31HG", y = "Expression")#+#tag="ENSG00000000003"
    p+facet_grid(rows = vars(box$项目名称), margins = TRUE)+ stat_compare_means(method = "kruskal.test", label.y =2 )
    # geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.2,fill = c("#999999"))
    #  Add p-value
    # aes(group="标本采集时间"),
    # p + stat_compare_means(method = "kruskal.test", label.y = 30 )
    # p + stat_compare_means(aes(group=cnv3),method = "kruskal.test", label.y = 500,comparisons =combn(1:3, 2, FUN = list) )+ #备注3
    #   stat_compare_means(comparisons = list(c("No changes","Shallow Deletion"),c("No changes","Deep Deletion"),c("Shallow Deletion","Deep Deletion")),method = "wilcox.test",label.y = as.numeric(seq(100,850,by=20)))+ #备注3
    #   stat_summary(aes(x = cnv3, y = number1),fun= "mean", geom = "point",position = position_dodge(0.75), shape = 23, size = 1, fill = "pink")
    # 
  }
  dev.off()
}


#1447 找问题
qqqqqqqqqqqq=resultALL_ref_sva[is.na(resultALL_ref_sva$项目结果),]
table(qqqqqqqqqqqq$项目名称)
qqqqqqqqqqqq=resultALL_remain[is.na(resultALL_remain$项目结果),]
table(qqqqqqqqqqqq$项目名称)
qqqqqqqqqqqq=resultALL_xxb[is.na(resultALL_xxb$项目结果),]
table(qqqqqqqqqqqq$项目名称)
qqqqqqqqqqqq=resultALL[is.na(resultALL$项目结果),]
table(qqqqqqqqqqqq$项目名称)
lapply(wuhu3,function(x)grepl("李满仓",unique(x$患者姓名)))
