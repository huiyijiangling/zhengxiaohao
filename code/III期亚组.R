rm(list=ls())
options(stringsAsFactors = F)
gc()
library(dplyr)
library(readxl)
library(tableone)
# JDW=read_excel("1421全胃13号最终版(TSN修改2).xlsx",sheet = 1,skip = 2)
# load("JDW_replustime.Rdata")
load("JDW_replustime(1).Rdata")
JDW=JDW_replustime
colnames(JDW)=stringr::str_split(colnames(JDW),"_",simplify = T)[,1]
JDW$rysj=stringr::str_split(JDW$rysj," ",simplify = T)[,1]
JDW$rysj1=JDW$rysj
dput(names(JDW))
JDW=subset(JDW,wwww==1)
JDW=subset(JDW,!is.na(JDW$pM))
str(JDW)
# load("JDW_replustime.Rdata")
# boxplot(JDW[,c("Age","sqts","shts")])
# 开始转变数据类型和清洗数据——方便做表
JDW$Sex=as.character(JDW$Sex)
JDW$Age=as.integer(JDW$Age)
JDW$Age2 = cut(JDW$Age,breaks = c(-Inf,65,Inf),right =T,ordered_result = T,
               labels = c("<65","≥65"))
JDW$zyts=as.integer(JDW$zyts)
# rysj、TS、Lauren、Borrmann、Diff列需要因子型排布
JDW$rysj=stringr::str_split(JDW$rysj,"-",simplify = T)[,1]
JDW$rysj=factor(JDW$rysj,
                levels = c("2014","2015","2016","2017",
                           "2018","2019","2020","2021"),
                labels = c("2014-2015","2014-2015","2016-2017",
                           "2016-2017","2018-2019","2018-2019",
                           "2020-2021", "2020-2021"),ordered = T)
JDW$rysj2=factor(JDW$rysj,
                 levels = c("2014-2015","2014-2015","2016-2017",
                            "2016-2017","2018-2019","2018-2019",
                            "2020-2021", "2020-2021"),
                 labels = c("Before2020","Before2020","Before2020",
                            "Before2020","Before2020","Before2020",
                            "After2020", "After2020"),ordered = T)

table(JDW$rysj)
JDW$neoadjuvant=factor(JDW$neoadjuvant,
                       levels = c("1","0"),
                       labels = c("1","0"),ordered = T)
JDW$TS=ifelse(is.na(JDW$TS),NA,
              ifelse(JDW$TS%in%c("0","1","2"),"proximal",
                     ifelse(JDW$TS%in%c("3","4","5"),"distal/total",NA)))
JDW$Lauren= ifelse(is.na(JDW$Lauren),NA,
                   ifelse(JDW$Lauren=="0","0",
                          ifelse(JDW$Lauren=="1","1",
                                 ifelse(JDW$Lauren=="2","2",NA))))
JDW$Lauren=factor(JDW$Lauren,levels = c(0:2),ordered = F,
                  labels=c("intestinal","mixed","diffuse"))
JDW$Borrmann=ifelse(is.na(JDW$Borrmann),NA,
                    ifelse(JDW$Borrmann%in%c("0","1"),"superficial",
                           ifelse(JDW$Borrmann%in%c("2","3","4"),"ulcerative",NA)))
JDW$Diff=factor(JDW$Diff,
                levels = c(0:5),
                labels = c("Poorly differentiated","Poorly differentiated",
                           "Poorly differentiated","Well differentiated",
                           "Well differentiated","Well differentiated"),ordered = T)

# 接着转变数据类型和清洗数据
JDW$lyvF=as.character(JDW$lyvF)
JDW$NF=as.character(JDW$NF)
JDW$signet=as.character(JDW$signet)
JDW$lyvF=ifelse(is.na(JDW$lyvF),NA,
                ifelse(JDW$lyvF=="0","Negative",
                       ifelse(JDW$lyvF=="1","Positive",NA)))
JDW$NF=ifelse(is.na(JDW$NF),NA,
              ifelse(JDW$NF=="0","Negative",
                     ifelse(JDW$NF=="1","Positive",NA)))
JDW$signet=factor(JDW$signet,levels = c(0:2),ordered = F,
                  labels = c("No-Signet","Partial-Signet","All-Signet"))
#AJCC 8版
JDW$pT1=ifelse(is.na(JDW$pT),NA,
               ifelse(JDW$pT=="0","0",
                      ifelse(JDW$pT%in%c("1","1a","1b"),"1",
                             ifelse(JDW$pT%in%c("2","2a","2b"),"2",
                                    ifelse(JDW$pT%in%c("3","3a","3b"),"3",
                                           ifelse(JDW$pT%in%c("4","4a"),"4",
                                                  ifelse(JDW$pT=="4b","5",NA)))))))
JDW$pN1=ifelse(is.na(JDW$pN),NA,
               ifelse(JDW$pN=="0","0",
                      ifelse(JDW$pN%in%c("1","1a","1b"),"1",
                             ifelse(JDW$pN%in%c("2","2a","2b"),"2",
                                    ifelse(JDW$pN%in%c("3","3a"),"3",
                                           ifelse(JDW$pN =="3b","5",NA))))))

JDW$tplusn=as.character(as.numeric(JDW$pT1)+as.numeric(JDW$pN1))
# JDW$pstage3 = ifelse(JDW$tplusn==6 & JDW$pT==4 , "IIIA",NA)

# 精细分期：直接手术病理分期和新辅助后的病理分期
JDW$pstage = ifelse(
  JDW$pM == "1","IV",
  ifelse(JDW$neoadjuvant=="0"&JDW$tplusn =="1", "IA",
         ifelse(JDW$neoadjuvant=="0"&JDW$tplusn =="2", "IB",
                ifelse(JDW$neoadjuvant=="0"&JDW$tplusn =="3", "IIA",
                       ifelse(JDW$neoadjuvant=="0"&JDW$tplusn =="4", "IIB",
                              ifelse(JDW$neoadjuvant=="0"&JDW$tplusn =="5", "IIIA",
                                     ifelse(JDW$neoadjuvant=="0"&JDW$tplusn=="6" & JDW$pT1=="4" , "IIIA",
                                            ifelse(JDW$neoadjuvant=="0"&JDW$tplusn%in%c("6","7"), "IIIB", 
                                                   ifelse(JDW$neoadjuvant=="0"&JDW$tplusn%in%c("8","9","10"), "IIIC",
                                                          ifelse(JDW$neoadjuvant=="1"&JDW$tplusn %in%c("1","2"), "I",
                                                                 ifelse(JDW$neoadjuvant=="1"&JDW$tplusn %in%c("3","4"), "II",
                                                                        ifelse(JDW$neoadjuvant=="1"&JDW$tplusn %in%c("5","6","7","8","9","10"),"III",
                                                                               ifelse(JDW$neoadjuvant=="1"&JDW$pT1=="1"&JDW$pN1=="5", "II",NA)))))))))))))

table(JDW$tplusn)
table(JDW$pstage)
#粗略分期
JDW$pstage1=ifelse(is.na(JDW$pstage),NA,
                   ifelse(JDW$pstage %in%c("IA", "IB","I"),"I",
                          ifelse(JDW$pstage %in%c("IIA","IIB","II"), "II",
                                 ifelse(JDW$pstage %in%c("IIIA","IIIB","IIIC","III"), "III",
                                        ifelse(JDW$pstage =="IV", "IV", NA)))))

JDW$pT=ifelse(is.na(JDW$pT),NA,
              ifelse(JDW$pT%in%c("0","1","1a","1b"),"early",
                     ifelse(JDW$pT%in%c("2","2a","2b",
                                        "3","3a","3b","4","4a","4b"),"advanced",NA)))
JDW$pN=ifelse(is.na(JDW$pN),NA,
              ifelse(JDW$pN=="0","Negative",
                     ifelse(JDW$pN%in%c("1","1a","1b","2",
                                        "2a","2b","3","3a","3b"),"Positive",NA)))
JDW$pM=ifelse(is.na(JDW$pM),NA,
              ifelse(JDW$pM=="0","Negative",
                     ifelse(JDW$pM%in%c("1","1a"),"Positive",NA)))
JDW$Smoking=as.character(JDW$Smoking)
JDW$Drinking=as.character(JDW$Drinking)
JDW$HD=as.character(JDW$HD)
JDW$hypertension=as.character(JDW$hypertension)
JDW$diabetes=as.character(JDW$diabetes)
JDW$Hyperlipidemia=as.character(JDW$Hyperlipidemia)
JDW$Antihypertensive=as.character(JDW$Antihypertensive)
JDW$antidiabetic=as.character(JDW$antidiabetic)
JDW$Statins=as.character(JDW$Statins)
JDW$Disease=ifelse((JDW$HD=="1")|(JDW$hypertension=="1")|(JDW$diabetes=="1")|(JDW$Hyperlipidemia=="1"),"Positive",
                   ifelse((JDW$HD=="0")|(JDW$hypertension=="0")|(JDW$diabetes=="0")|(JDW$Hyperlipidemia=="0"),"Negative","NA"))
table(JDW$Disease)
JDW$Tumormargin=ifelse(is.na(JDW$Tumormargin),NA,
                       ifelse(JDW$Tumormargin=="R0","Negative",
                              ifelse(JDW$Tumormargin%in%c("R1","R2"),"Positive","Wrong")))
JDW$surgerymargin=ifelse(is.na(JDW$surgerymargin),NA,
                         ifelse(JDW$surgerymargin=="R0","Negative",
                                ifelse(JDW$surgerymargin%in%c("R1","R2"),"Positive","Wrong")))
JDW$sqts=as.integer(JDW$sqts)
JDW$shts=as.integer(JDW$shts)
JDW$yqts=as.integer(JDW$yqts)
JDW$qzryts=as.integer(JDW$qzryts)
JDW$qzssts=as.integer(JDW$qzssts)
JDW$dsts=as.integer(JDW$dsts)
JDW$NAC=as.character(JDW$NAC)
JDW$NAR=as.character(JDW$NAR)
JDW$PAC=as.character(JDW$PAC)
JDW$PAR=as.character(JDW$PAR)
dput(names(JDW)) # 输出据集变量名称
str(JDW)
quantile(JDW$dsts)
JDW$dsts2=cut(JDW$dsts,breaks = c(-Inf,14,30,Inf),right =T,ordered_result = T,
              labels = c("WT≤14","14<WT≤30","30<WT"))


box=JDW[JDW$pstage1%in%c("III","IV"),]#晚期胃癌数据集

tabledsts <- CreateTableOne(vars =  c("Sex", "Age","Age2","dsts2","Smoking", "Drinking","Disease","neoadjuvant","TS",  "size", "Lauren", "Borrmann",
                                      "Diff", "lyvF", "NF","signet", "pT", "pN", "pM", "pstage1","Tumormargin","surgerymargin"),
                            data = box,includeNA = F, addOverall = T,strata = "xjdsts3")
tabledsts1=print(tabledsts,  formatOptions = list(big.mark = ","),
                 showAllLevels = TRUE, quote =  F, noSpaces = T,
                 missing=T, printToggle = T,explain=F)
write.csv(tabledsts1,quote=T,file="III期亚组/III期TTS切值30.csv")



