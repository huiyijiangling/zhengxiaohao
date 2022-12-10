library(dplyr)
library(readxl)
library(tableone)
JDW_replustime$hypertension
Sex, Age, TS, size, Lauren, 
  Borrmann, Diff, lyvF, NF, signet, pT, pN, pM
  Smoking, Drinking, HD, hypertension, diabetes
box_completecase$Sex=as.character(box_completecase$Sex)
box_completecase$Age=as.integer(box_completecase$Age)
box_completecase$Age2 = cut(box_completecase$Age,breaks = c(-Inf,65,Inf),right =T,ordered_result = T,
               labels = c("<65","≥65"))
box_completecase$zyts=as.integer(box_completecase$zyts)
# rysj、TS、Lauren、Borrmann、Diff列需要因子型排布
box_completecase$rysj=stringr::str_split(box_completecase$rysj,"-",simplify = T)[,1]
box_completecase$rysj=factor(box_completecase$rysj,
                levels = c("2014","2015","2016","2017",
                           "2018","2019","2020","2021"),
                labels = c("2014-2015","2014-2015","2016-2017",
                           "2016-2017","2018-2019","2018-2019",
                           "2020-2021", "2020-2021"),ordered = T)
table(box_completecase$rysj)
box_completecase$COVID=ifelse(is.na(box_completecase$rysj1),NA,
                 ifelse(box_completecase$rysj1<"2020-01-22","COVID-19B",
                        ifelse(box_completecase$rysj1>="2020-01-22","COVID-19A",box_completecase$rysj1)))
box_completecase$COVID=factor(box_completecase$COVID,
                 levels = c("COVID-19B","COVID-19A"),
                 labels = c("COVID-19B","COVID-19A"),ordered = T)
box_completecase$neoadjuvant=factor(box_completecase$neoadjuvant,
                       levels = c("1","0"),
                       labels = c("有","无"),ordered = T)
box_completecase$TS=ifelse(is.na(box_completecase$TS),NA,
              ifelse(box_completecase$TS%in%c("0","1","2"),"proximal",
                     ifelse(box_completecase$TS%in%c("3","4","5"),"distal/total",NA)))
box_completecase$Lauren= ifelse(is.na(box_completecase$Lauren),NA,
                   ifelse(box_completecase$Lauren=="0","0",
                          ifelse(box_completecase$Lauren=="1","1",
                                 ifelse(box_completecase$Lauren=="2","2",NA))))
box_completecase$Lauren=factor(box_completecase$Lauren,levels = c(0:2),ordered = F,
                  labels=c("intestinal","mixed","diffuse"))
box_completecase$Borrmann=ifelse(is.na(box_completecase$Borrmann),NA,
                    ifelse(box_completecase$Borrmann%in%c("0","1"),"superficial",
                           ifelse(box_completecase$Borrmann%in%c("2","3","4"),"ulcerative",NA)))
box_completecase$Diff=factor(box_completecase$Diff,
                levels = c(0:5),
                labels = c("Poorly differentiated","Poorly differentiated",
                           "Poorly differentiated","Well differentiated",
                           "Well differentiated","Well differentiated"),ordered = T)


# 接着转变数据类型和清洗数据
box_completecase$lyvF=as.character(box_completecase$lyvF)
box_completecase$NF=as.character(box_completecase$NF)
box_completecase$signet=as.character(box_completecase$signet)
box_completecase$lyvF=ifelse(is.na(box_completecase$lyvF),NA,
                ifelse(box_completecase$lyvF=="0","Negative",
                       ifelse(box_completecase$lyvF=="1","Positive",NA)))
box_completecase$NF=ifelse(is.na(box_completecase$NF),NA,
              ifelse(box_completecase$NF=="0","Negative",
                     ifelse(box_completecase$NF=="1","Positive",NA)))
box_completecase$signet=factor(box_completecase$signet,levels = c(0:2),ordered = F,
                  labels = c("No-Signet","Partial-Signet","All-Signet"))
#AJCC 8版
box_completecase$pT1=ifelse(is.na(box_completecase$pT),NA,
               ifelse(box_completecase$pT=="0","0",
                      ifelse(box_completecase$pT%in%c("1","1a","1b"),"1",
                             ifelse(box_completecase$pT%in%c("2","2a","2b"),"2",
                                    ifelse(box_completecase$pT%in%c("3","3a","3b"),"3",
                                           ifelse(box_completecase$pT%in%c("4","4a"),"4",
                                                  ifelse(box_completecase$pT=="4b","5",NA)))))))
box_completecase$pN1=ifelse(is.na(box_completecase$pN),NA,
               ifelse(box_completecase$pN=="0","0",
                      ifelse(box_completecase$pN%in%c("1","1a","1b"),"1",
                             ifelse(box_completecase$pN%in%c("2","2a","2b"),"2",
                                    ifelse(box_completecase$pN%in%c("3","3a"),"3",
                                           ifelse(box_completecase$pN =="3b","5",NA))))))

box_completecase$tplusn=as.character(as.numeric(box_completecase$pT1)+as.numeric(box_completecase$pN1))
# box_completecase$pstage3 = ifelse(box_completecase$tplusn==6 & box_completecase$pT==4 , "IIIA",NA)
box_completecase$pstage = ifelse(
  box_completecase$pM == "1","IV",
  ifelse(box_completecase$tplusn =="1", "IA",
         ifelse(box_completecase$tplusn =="2", "IB",
                ifelse(box_completecase$tplusn =="3", "IIA",
                       ifelse(box_completecase$tplusn =="4", "IIB",
                              ifelse(box_completecase$tplusn =="5", "IIIA",
                                     ifelse(box_completecase$tplusn=="6" & box_completecase$pT1=="4" , "IIIA",
                                            ifelse(box_completecase$tplusn%in%c("6","7"), "IIIB", 
                                                   ifelse(box_completecase$tplusn%in%c("8","9","10"), "IIIC", NA)))))))))
table(box_completecase$tplusn)
table(box_completecase$pstage)
box_completecase$pstage1=ifelse(is.na(box_completecase$pstage),NA,
                   ifelse(box_completecase$pstage %in%c("IA", "IB"),"I",
                          ifelse(box_completecase$pstage %in%c("IIA","IIB"), "II",
                                 ifelse(box_completecase$pstage %in%c("IIIA","IIIB","IIIC"), "III",
                                        ifelse(box_completecase$pstage =="IV", "IV",NA)))))


box_completecase$pT=ifelse(is.na(box_completecase$pT),NA,
              ifelse(box_completecase$pT%in%c("0","1","1a","1b"),"early",
                     ifelse(box_completecase$pT%in%c("2","2a","2b",
                                        "3","3a","3b","4","4a","4b"),"advanced",NA)))
box_completecase$pN=ifelse(is.na(box_completecase$pN),NA,
              ifelse(box_completecase$pN=="0","Negative",
                     ifelse(box_completecase$pN%in%c("1","1a","1b","2",
                                        "2a","2b","3","3a","3b"),"Positive",NA)))
box_completecase$pM=ifelse(is.na(box_completecase$pM),NA,
              ifelse(box_completecase$pM=="0","Negative",
                     ifelse(box_completecase$pM%in%c("1","1a"),"Positive",NA)))
box_completecase$Smoking=as.character(box_completecase$Smoking)
box_completecase$Drinking=as.character(box_completecase$Drinking)
box_completecase$HD=as.character(box_completecase$HD)
box_completecase$hypertension=as.character(box_completecase$hypertension)
box_completecase$diabetes=as.character(box_completecase$diabetes)
box_completecase$Hyperlipidemia=as.character(box_completecase$Hyperlipidemia)
box_completecase$Antihypertensive=as.character(box_completecase$Antihypertensive)
box_completecase$antidiabetic=as.character(box_completecase$antidiabetic)
box_completecase$Statins=as.character(box_completecase$Statins)
box_completecase$Disease=ifelse((box_completecase$HD=="1")|(box_completecase$hypertension=="1")|(box_completecase$diabetes=="1")|(box_completecase$Hyperlipidemia=="1"),"Positive",
                   ifelse((box_completecase$HD=="0")|(box_completecase$hypertension=="0")|(box_completecase$diabetes=="0")|(box_completecase$Hyperlipidemia=="0"),"Negative","NA"))
table(box_completecase$Disease)
box_completecase$Tumormargin=ifelse(is.na(box_completecase$Tumormargin),NA,
                       ifelse(box_completecase$Tumormargin=="R0","Negative",
                              ifelse(box_completecase$Tumormargin%in%c("R1","R2"),"Positive","Wrong")))
box_completecase$surgerymargin=ifelse(is.na(box_completecase$surgerymargin),NA,
                         ifelse(box_completecase$surgerymargin=="R0","Negative",
                                ifelse(box_completecase$surgerymargin%in%c("R1","R2"),"Positive","Wrong")))