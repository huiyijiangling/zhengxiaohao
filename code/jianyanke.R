#strptime,strptime,as.POSIXct,as.POSIXlt本质都不好
#strftime 最好
#NULL 为空
#先找 空 na null
#table(grep("NULL",resultALL[[5]]$报告时间)) 报告时间可以是空
#报告格式不统一找na
# deparse(substitute(box))
# eval(parse(text = "box"))
library(data.table)
library(doParallel)
library(parallel)#
library(snow)#
library(dplyr)
library(lubridate)
library(readxl)
library(tidyr)
library(sva)
library(ggplot2)
library(ggpubr)
library(Publish)
library(nortest)
rm(list=ls())
options(stringsAsFactors = F)
gc()

filenames_1 <- list.files('./jianyan/',pattern = ".txt",full.name=TRUE)

resultALL1 <- lapply(filenames_1[1], function(fl)   data.table::fread(file =fl, header = T, na.strings = "NULL",sep = "\t",verbose = F, integer64 = 'numeric',encoding =  "unknown",fill = F))
resultALL2 <- lapply(filenames_1[2:8], function(fl)   data.table::fread(file =fl, header = T, na.strings = "NULL",sep = "\t",verbose = F, integer64 = 'numeric',encoding =  "UTF-8",fill = T))

check_chuanhang=lapply(1:7,function(x) resultALL2[[x]][is.na(resultALL2[[x]]$就诊卡号),]$申请科室)
# kkkk=paste(resultALL2[[6]][grep("渗出液",resultALL2[[6]]$申请科室)-1,]$参考值,resultALL2[[6]][grep("渗出液",resultALL2[[6]]$申请科室),]$申请科室,sep=" ")
# resultALL2=lapply(1:7,function(x) 
lengths(check_chuanhang)
# 0   0   0   0 168 168 114
for (x in 5:7) {
  resultALL2[[x]][grep("渗出液",resultALL2[[x]]$申请科室)-1,15] <-paste(resultALL2[[x]][grep("渗出液",resultALL2[[x]]$申请科室)-1,]$参考值,resultALL2[[x]][grep("渗出液",resultALL2[[x]]$申请科室),]$申请科室,sep=" ")
  resultALL2[[x]]=subset(resultALL2[[x]],!grepl("渗出液",resultALL2[[x]]$申请科室))
}   
resultALL=c(resultALL1,resultALL2)
resultALL=lapply(resultALL,as.data.frame)
# resultALL=data.table::rbindlist(resultALL,use.names = T)#报错             
class(resultALL[[1]])                 
Encoding(colnames(resultALL[[1]]))


cl <- makeCluster(8)
registerDoParallel(cl)
unknowntoutf8 <- function(x){
  x=resultALL[[x]]
  a=enc2utf8(rownames(x))
  b=enc2utf8(colnames(x))
  x=snow::parApply(cl,x,2,function(x) enc2utf8(as.character(x)))
  x=as.data.frame(x)
  rownames(x)=a
  colnames(x)=b
  return(x)
}
resultALLss=list()
for(i in 1:8){
  resultALLss[[i]]=unknowntoutf8(i)
}
stopCluster(cl)

resultALL=resultALLss
resultALL=lapply(resultALL,as.data.frame)
Encoding(resultALL[[5]]$标本采集时间)
Encoding(resultALL[[5]]$报告时间)
Encoding(resultALL1[[1]]$标本采集时间)
Encoding(resultALL1[[1]]$报告时间)
table(is.na(resultALL[[5]]$标本采集时间))
resultALL[[1]]$标本采集时间=gsub("\\/","-",resultALL[[1]]$标本采集时间)
resultALL[[1]]$报告时间=gsub("\\/","-",resultALL[[1]]$报告时间)
resultALL[[1]]$报告时间 <- lubridate::ymd_hm(resultALL[[1]]$报告时间,tz=Sys.timezone(location = TRUE))
resultALL[[1]]$标本采集时间 <- lubridate::ymd_hm(resultALL[[1]]$标本采集时间,tz=Sys.timezone(location = TRUE))
for (x in 2:8) {
  resultALL[[x]]$报告时间 <- lubridate::ymd_hms(resultALL[[x]]$报告时间,tz=Sys.timezone(location = TRUE))
  resultALL[[x]]$标本采集时间 <- lubridate::ymd_hms(resultALL[[x]]$标本采集时间,tz=Sys.timezone(location = TRUE))
}
resultALL=data.table::rbindlist(resultALL,use.names = T,fill = F)

table(is.na(resultALL$报告时间))
table(is.na(resultALL$标本采集时间))
table(grep("NULL",resultALL$报告时间))
table(grep("NULL",resultALL$标本采集时间))
table(grep("NULL",resultALL$就诊卡号))
table(grepl("NULL",resultALL$病案号))

resultALL_colnames=colnames(resultALL)
cl <- makeCluster(15)
registerDoParallel(cl)
resultALL=foreach(x=resultALL)%dopar%gsub("NULL",NA,x)
# parApply(cl,resultALL,2,function(x)gsub("NULL",NA,x))#其实在na.string里面有问题
stopCluster(cl)
resultALL=as.data.frame(resultALL,col.names = resultALL_colnames)

table(is.na(resultALL$就诊卡号))
class(resultALL$就诊卡号)
resultALL$病案号=as.numeric(resultALL$病案号)
resultALL$就诊卡号=as.numeric(resultALL$就诊卡号)
table(is.na(resultALL$病案号))
table(is.na(resultALL$就诊卡号))

resultALL$报告时间 <- lubridate::ymd_hms(resultALL$报告时间,tz=Sys.timezone(location = TRUE))
resultALL$标本采集时间 <- lubridate::ymd_hms(resultALL$标本采集时间,tz=Sys.timezone(location = TRUE))

save(resultALL,file="resultALL01.Rdata")
# load(file="resultALL01.Rdata")

if(T){
wuhuqifei=resultALL
wuhuqifei$IDID=seq(1,nrow(wuhuqifei),by=1)
wuhuqifei$项目结果=as.numeric(wuhuqifei$项目结果)
wuhuqifei=subset(wuhuqifei,is.na(wuhuqifei$项目结果))
wuhuqifei=resultALL[wuhuqifei$IDID,]
# manual check
huluwa=as.data.frame(table(wuhuqifei$项目结果,wuhuqifei$项目名称,useNA = "ifany"))
huluwa=subset(huluwa,huluwa$Freq!=0)  
write.csv(huluwa,file = "huluwa_check.csv")
wuhu_reduced=readxl::read_xlsx(path =r"(C:\Users\zxh\Desktop\R\jianyanke\huluwa_check.xlsx)",sheet =1 )
wuhuqifei=subset(wuhuqifei,wuhuqifei$项目名称%in%wuhu_reduced$Var2)
unique(wuhu_reduced$Var2)
unique(wuhuqifei$项目名称)
#ok
wuhuqifei=subset(wuhuqifei,select = c("就诊卡号","标本采集时间"))
wuhuqifei=distinct(wuhuqifei)
wuhuqifei=unite(wuhuqifei,col="ID_time","就诊卡号","标本采集时间",sep = "_",remove = F)
dput(names(resultALL))
# c("申请科室", "患者姓名", "患者性别", "患者年龄", "就诊卡号", 
#   "病案号", "诊断", "标本采集时间", "专业名称", "报告时间", "项目名称", 
#   "项目中文注释", "项目结果", "单位", "参考值")
resultALL=unite(resultALL,col="ID_time","就诊卡号","标本采集时间",sep = "_",remove = F)
dput(names(resultALL))
# c("申请科室", "患者姓名", "患者性别", "患者年龄", "ID_time",
#   "就诊卡号", "病案号", "诊断", "标本采集时间", "专业名称", "报告时间",
#   "项目名称", "项目中文注释", "项目结果", "单位", "参考值")
wuhu2=subset(resultALL,ID_time %in%wuhuqifei$ID_time)
resultALL_no_na=subset(resultALL,!ID_time %in%wuhuqifei$ID_time)
# wuhu2=merge(resultALL,wuhuqifei,by = c("就诊卡号","标本采集时间"))
# wuhu2=unite(wuhu2,col="ID_time","就诊卡号","标本采集时间",sep = "_",remove = F)
wuhu3=split(wuhu2,wuhu2$ID_time)
if(F){
names(wuhu3)=NULL
for (i in 1:length(wuhu3)){
openxlsx::write.xlsx(wuhu3[[i]],file=paste0("./outlier test/",i,".xlsx"))}
}
wuhu3_after_check <- list.files('./outlier test version 20210130/',pattern = ".xlsx",full.name=TRUE)
wuhu3_after_check=gtools::mixedsort(wuhu3_after_check)
wuhu3_after_check <- unlist(lapply(wuhu3_after_check, function(fl) readxl::read_xlsx(fl)$ID_time[1]))
duplicated(wuhu3_after_check)
duplicated(wuhu3)
wuhu3=wuhu3[match(wuhu3_after_check,names(wuhu3))]
wuhu3=foreach(i=1:length(wuhu3))%do%{
  wuhu3[[i]]$项目结果=as.numeric(wuhu3[[i]]$项目结果)
  return(wuhu3[[i]])
}
}
# match 是用 匹配,需要单独的，前方的是你想匹配的东西（想要的顺序），后方是原顺序，待调整的放在后方
#导入已经处理好的 101 个手动检查的结果 对齐结果 防止错误
lapply(1:length(wuhu3),function(i) wuhu3[[i]][is.na(wuhu3[[i]]$项目结果),]$项目名称)
# 手动处理异常值 从逻辑上是后面先算好再提回来算的
# wuhu3[[1]][which(wuhu3[[1]]$项目名称=="NRBC%"),"项目结果"]=NA
if(T){
  #并不一定真实值 只是按照上限值来插补
for(i in 1:12){
  wuhu3[[i]][which(wuhu3[[i]]$项目名称=="NRBC%"),"项目结果"]=1.9
}
  
  wuhu3[[13]][which(wuhu3[[13]]$项目名称=="AT-III"),"项目结果"]=2.0
  
  wuhu3[[14]][which(wuhu3[[14]]$项目名称=="TT"),"项目结果"]=250

  wuhu3[[15]][which(wuhu3[[15]]$项目名称=="TT"),"项目结果"]=250
  wuhu3[[15]][which(wuhu3[[15]]$项目名称=="PT(s)"),"项目结果"]=450
  wuhu3[[15]][which(wuhu3[[15]]$项目名称=="PT(INR)"),"项目结果"]=39.99
  wuhu3[[15]][which(wuhu3[[15]]$项目名称=="PT(a)"),"项目结果"]=1.5
  wuhu3[[15]][which(wuhu3[[15]]$项目名称=="PT(r)"),"项目结果"]=40
  wuhu3[[15]][which(wuhu3[[15]]$项目名称=="D-D"),"项目结果"]=250
  wuhu3[[15]][which(wuhu3[[15]]$项目名称=="APTT"),"项目结果"]=50.0

  wuhu3[[16]][which(wuhu3[[16]]$项目名称=="PT(s)"),"项目结果"]=450
  wuhu3[[16]][which(wuhu3[[16]]$项目名称=="PT(a)"),"项目结果"]=1.5
  wuhu3[[16]][which(wuhu3[[16]]$项目名称=="PT(r)"),"项目结果"]=40
  wuhu3[[16]][which(wuhu3[[16]]$项目名称=="PT(INR)"),"项目结果"]=39.99
  wuhu3[[16]][which(wuhu3[[16]]$项目名称=="TT"),"项目结果"]=250
  wuhu3[[16]][which(wuhu3[[16]]$项目名称=="FIB"),"项目结果"]=0.15
  wuhu3[[16]][which(wuhu3[[16]]$项目名称=="APTT"),"项目结果"]=200

  wuhu3[[17]][which(wuhu3[[17]]$项目名称=="GA%"),"项目结果"]=14.3#median随机
  
  wuhu3[[18]][which(wuhu3[[18]]$项目名称=="TT"),"项目结果"]=250

  wuhu3[[19]][which(wuhu3[[19]]$项目名称=="TT"),"项目结果"]=250
  wuhu3[[19]][which(wuhu3[[19]]$项目名称=="APTT"),"项目结果"]=200

  wuhu3[[20]][which(wuhu3[[20]]$项目名称=="TT"),"项目结果"]=13.4
  
  wuhu3[[21]][which(wuhu3[[21]]$项目名称=="TT"),"项目结果"]=11.5

  wuhu3[[22]][which(wuhu3[[22]]$项目名称=="TT"),"项目结果"]=13.5

  wuhu3[[23]][which(wuhu3[[23]]$项目名称=="TT"),"项目结果"]=26.6

  wuhu3[[24]][which(wuhu3[[24]]$项目名称=="TT"),"项目结果"]=12.7

  wuhu3[[25]][which(wuhu3[[25]]$项目名称=="TT"),"项目结果"]=14.0

  wuhu3[[26]][which(wuhu3[[26]]$项目名称=="PALB"),"项目结果"]=23.0
  wuhu3[[26]][which(wuhu3[[26]]$项目名称=="LPa"),"项目结果"]=80.0
  wuhu3[[26]][which(wuhu3[[26]]$项目名称=="LDL-CHO"),"项目结果"]=4.0
  wuhu3[[26]][which(wuhu3[[26]]$项目名称=="DBIL"),"项目结果"]=0.1
  wuhu3[[26]][nrow(wuhu3[[26]])+1,]=wuhu3[[26]][which(wuhu3[[26]]$项目名称=="DBIL"),]
  wuhu3[[26]][nrow(wuhu3[[26]]),c("项目名称","项目中文注释","项目结果","参考值")]=c("IBIL","间接(游离)胆红素",0.2,"0-11.97")
  ########
  
  #错的？？？
  wuhu3[[27]][which(wuhu3[[27]]$项目名称=="APTT"),"项目结果"]=27.5
  wuhu3[[27]][which(wuhu3[[27]]$项目名称=="TT"),"项目结果"]=17
  wuhu3[[27]][which(wuhu3[[27]]$项目名称=="D-D"),"项目结果"]=0.6
  wuhu3[[27]][which(wuhu3[[27]]$项目名称=="FIB"),"项目结果"]=3.3
  wuhu3[[27]][which(wuhu3[[27]]$项目名称=="AT-III"),"项目结果"]=92.6

  #
  for(i in 28:34){
    wuhu3[[i]][which(wuhu3[[i]]$项目名称=="APTT"),"项目结果"]=250
  }

  wuhu3[[35]][which(wuhu3[[35]]$项目名称=="D-D"),"项目结果"]=500

  wuhu3[[36]][which(wuhu3[[36]]$项目名称=="TT"),"项目结果"]=200

  wuhu3[[37]][which(wuhu3[[37]]$项目名称=="APTT"),"项目结果"]=250

  wuhu3[[38]][which(wuhu3[[38]]$项目名称=="APTT"),"项目结果"]=250
  
  wuhu3[[39]][which(wuhu3[[39]]$项目名称=="TT"),"项目结果"]=200

  wuhu3[[40]][which(wuhu3[[40]]$项目名称=="NEUT%"),"项目结果"]=1
  wuhu3[[40]][which(wuhu3[[40]]$项目名称=="NEUT#"),"项目结果"]=0.004
  wuhu3[[40]][which(wuhu3[[40]]$项目名称=="EO%"),"项目结果"]=31.5
  wuhu3[[40]][which(wuhu3[[40]]$项目名称=="EO#"),"项目结果"]=0.126

  wuhu3[[41]][which(wuhu3[[41]]$项目名称=="D-D"),"项目结果"]=500
  wuhu3[[41]][which(wuhu3[[41]]$项目名称=="FDP"),"项目结果"]=15000

  wuhu3[[42]][which(wuhu3[[42]]$项目名称=="APTT"),"项目结果"]=250
  wuhu3[[42]][which(wuhu3[[42]]$项目名称=="TT"),"项目结果"]=200
  
  wuhu3[[43]][which(wuhu3[[43]]$项目名称=="FDP"),"项目结果"]=300

  wuhu3[[44]][which(wuhu3[[44]]$项目名称=="TT"),"项目结果"]=19.1

  wuhu3[[45]][which(wuhu3[[45]]$项目名称=="APTT"),"项目结果"]=250

  wuhu3[[46]][which(wuhu3[[46]]$项目名称=="APTT"),"项目结果"]=250

  wuhu3[[47]][which(wuhu3[[47]]$项目名称=="P-LCR"),"项目结果"]=41.2
  wuhu3[[47]][which(wuhu3[[47]]$项目名称=="FDP"),"项目结果"]=300
  wuhu3[[47]][which(wuhu3[[47]]$项目名称=="D-D"),"项目结果"]=500

  wuhu3[[48]][which(wuhu3[[48]]$项目名称=="PT(a)"),"项目结果"]=3.5
  wuhu3[[48]][which(wuhu3[[48]]$项目名称=="PT(r)"),"项目结果"]=11.0
  wuhu3[[48]][which(wuhu3[[48]]$项目名称=="PT(INR)"),"项目结果"]=10.99
  wuhu3[[48]][which(wuhu3[[48]]$项目名称=="D-D"),"项目结果"]=500
  wuhu3[[48]][which(wuhu3[[48]]$项目名称=="TT"),"项目结果"]=200
  wuhu3[[48]][which(wuhu3[[48]]$项目名称=="APTT"),"项目结果"]=250
  wuhu3[[48]][which(wuhu3[[48]]$项目名称=="PT(s)"),"项目结果"]=250
  wuhu3[[48]][which(wuhu3[[48]]$项目名称=="FIB"),"项目结果"]=0.1
  
  wuhu3[[49]][which(wuhu3[[49]]$项目名称=="PT(r)"),"项目结果"]=11.0
  wuhu3[[49]][which(wuhu3[[49]]$项目名称=="PT(INR)"),"项目结果"]=10.99
  wuhu3[[49]][which(wuhu3[[49]]$项目名称=="PT(a)"),"项目结果"]=3.5
  wuhu3[[49]][which(wuhu3[[49]]$项目名称=="PT(s)"),"项目结果"]=250
  wuhu3[[49]][which(wuhu3[[49]]$项目名称=="FIB"),"项目结果"]=0.1
  
  wuhu3[[50]][which(wuhu3[[50]]$项目名称=="PT(a)"),"项目结果"]=3.5
  wuhu3[[50]][which(wuhu3[[50]]$项目名称=="PT(r)"),"项目结果"]=11.0
  wuhu3[[50]][which(wuhu3[[50]]$项目名称=="PT(INR)"),"项目结果"]=10.99
  wuhu3[[50]][which(wuhu3[[50]]$项目名称=="TT"),"项目结果"]=200
  wuhu3[[50]][which(wuhu3[[50]]$项目名称=="APTT"),"项目结果"]=250
  wuhu3[[50]][which(wuhu3[[50]]$项目名称=="PT(s)"),"项目结果"]=250
  wuhu3[[50]][which(wuhu3[[50]]$项目名称=="FIB"),"项目结果"]=0.1

  wuhu3[[51]][which(wuhu3[[51]]$项目名称=="FIB"),"项目结果"]=0.1
  wuhu3[[51]][which(wuhu3[[51]]$项目名称=="D-D"),"项目结果"]=500

  wuhu3[[52]][which(wuhu3[[52]]$项目名称=="D-D"),"项目结果"]=500

  wuhu3[[53]][which(wuhu3[[53]]$项目名称=="D-D"),"项目结果"]=120
  
  wuhu3[[54]][which(wuhu3[[54]]$项目名称=="APTT"),"项目结果"]=250

  wuhu3[[55]][which(wuhu3[[55]]$项目名称=="FIB"),"项目结果"]=3.3
  
  wuhu3[[56]][which(wuhu3[[56]]$项目名称=="FIB"),"项目结果"]=2.9

  wuhu3[[57]][which(wuhu3[[57]]$项目名称=="APTT"),"项目结果"]=250
  wuhu3[[58]][which(wuhu3[[58]]$项目名称=="APTT"),"项目结果"]=250
  wuhu3[[59]][which(wuhu3[[59]]$项目名称=="APTT"),"项目结果"]=250

  wuhu3[[60]][which(wuhu3[[60]]$项目名称=="PT(a)"),"项目结果"]=3.5#<2018-01-08
  wuhu3[[60]][which(wuhu3[[60]]$项目名称=="PT(r)"),"项目结果"]=11.0
  wuhu3[[60]][which(wuhu3[[60]]$项目名称=="PT(INR)"),"项目结果"]=10.99
  wuhu3[[60]][which(wuhu3[[60]]$项目名称=="PT(s)"),"项目结果"]=250

  wuhu3[[61]][which(wuhu3[[61]]$项目名称=="PT(a)"),"项目结果"]=3.5#<2018-01-08
  wuhu3[[61]][which(wuhu3[[61]]$项目名称=="PT(r)"),"项目结果"]=11.0
  wuhu3[[61]][which(wuhu3[[61]]$项目名称=="PT(INR)"),"项目结果"]=10.99
  wuhu3[[61]][which(wuhu3[[61]]$项目名称=="PT(s)"),"项目结果"]=250
  
  wuhu3[[62]][which(wuhu3[[62]]$项目名称=="TT"),"项目结果"]=15.5
  wuhu3[[62]][which(wuhu3[[62]]$项目名称=="APTT"),"项目结果"]=24.8

  # 2018-10-18 16:55:47
  wuhu3[[63]][which(wuhu3[[63]]$项目名称=="PT(a)"),"项目结果"]=1.5#>2018-01-08
  wuhu3[[63]][which(wuhu3[[63]]$项目名称=="PT(r)"),"项目结果"]=40#>2018-01-08
  wuhu3[[63]][which(wuhu3[[63]]$项目名称=="PT(INR)"),"项目结果"]=39.99#>2018-01-08
  wuhu3[[63]][which(wuhu3[[63]]$项目名称=="PT(s)"),"项目结果"]=250#<2019-6-24

  wuhu3[[64]][which(wuhu3[[64]]$项目名称=="HFLC%"),"项目结果"]=0
  wuhu3[[64]][which(wuhu3[[64]]$项目名称=="IG%"),"项目结果"]=1.3#>2018-01-08
  wuhu3[[64]][which(wuhu3[[64]]$项目名称=="NRBC%"),"项目结果"]=0.1#>2018-01-08
  
  wuhu3[[65]][which(wuhu3[[65]]$项目名称=="APTT"),"项目结果"]=250#<2019-06-24

  wuhu3[[66]][which(wuhu3[[66]]$项目名称=="NRBC%"),"项目结果"]=0
  
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="HFLC%"),"项目结果"]=0.0
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="IG%"),"项目结果"]=0.0
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="NEUT%"),"项目结果"]=0
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="NEUT#"),"项目结果"]=0
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="LYMPH%"),"项目结果"]=75
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="LYMPH#"),"项目结果"]=0.06
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="MONO%"),"项目结果"]=20
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="MONO#"),"项目结果"]=0.02
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="EO%"),"项目结果"]=0
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="EO#"),"项目结果"]=0
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="BASO%"),"项目结果"]=0
  wuhu3[[67]][which(wuhu3[[67]]$项目名称=="BASO#"),"项目结果"]=0

  wuhu3[[68]][which(wuhu3[[68]]$项目名称=="APTT"),"项目结果"]=62

  wuhu3[[69]][which(wuhu3[[69]]$项目名称=="TT"),"项目结果"]=200#<2019-06-24
  wuhu3[[69]][which(wuhu3[[69]]$项目名称=="APTT"),"项目结果"]=250#<2019-06-24

  wuhu3[[70]][which(wuhu3[[70]]$项目名称=="TT"),"项目结果"]=200#<2019-06-24

  wuhu3[[71]][which(wuhu3[[71]]$项目名称=="HFLC%"),"项目结果"]=0.1
  wuhu3[[71]][which(wuhu3[[71]]$项目名称=="IG%"),"项目结果"]=1.1#>2018-01-08
  wuhu3[[71]][which(wuhu3[[71]]$项目名称=="NRBC%"),"项目结果"]=0#>2018-01-08

  wuhu3[[72]][which(wuhu3[[72]]$项目名称=="GA%"),"项目结果"]=13.8#median
  
  wuhu3[[73]][which(wuhu3[[73]]$项目名称=="APTT"),"项目结果"]=250
  
  wuhu3[[74]][which(wuhu3[[74]]$项目名称=="APTT"),"项目结果"]=30

  wuhu3[[75]][which(wuhu3[[75]]$项目名称=="TT"),"项目结果"]=18.8
  wuhu3[[75]][which(wuhu3[[75]]$项目名称=="APTT"),"项目结果"]=22.1
  
  wuhu3[[76]][which(wuhu3[[76]]$项目名称=="HFLC%"),"项目结果"]=0
  wuhu3[[76]][which(wuhu3[[76]]$项目名称=="IG%"),"项目结果"]=0.4
  wuhu3[[76]][which(wuhu3[[76]]$项目名称=="P-LCR"),"项目结果"]=39.8
  wuhu3[[76]][which(wuhu3[[76]]$项目名称=="NRBC%"),"项目结果"]=0

  wuhu3[[77]][which(wuhu3[[77]]$项目名称=="APTT"),"项目结果"]=26.8

  wuhu3[[78]][which(wuhu3[[78]]$项目名称=="APTT"),"项目结果"]=250#<2019-06-24

  wuhu3[[79]][which(wuhu3[[79]]$项目名称=="TT"),"项目结果"]=200#<2019-06-24

  wuhu3[[80]][which(wuhu3[[80]]$项目名称=="APTT"),"项目结果"]=39

  wuhu3[[81]][which(wuhu3[[81]]$项目名称=="APTT"),"项目结果"]=250#<2019-06-24
  wuhu3[[81]][which(wuhu3[[81]]$项目名称=="PT(a)"),"项目结果"]=1.5#>2018-01-08
  wuhu3[[81]][which(wuhu3[[81]]$项目名称=="PT(r)"),"项目结果"]=40#>2018-01-08
  wuhu3[[81]][which(wuhu3[[81]]$项目名称=="PT(INR)"),"项目结果"]=39.99#>2018-01-08
  wuhu3[[81]][which(wuhu3[[81]]$项目名称=="PT(s)"),"项目结果"]=250#<2019-6-24
  
  wuhu3[[82]][which(wuhu3[[82]]$项目名称=="APTT"),"项目结果"]=250#<2019-06-24

  wuhu3[[83]][which(wuhu3[[83]]$项目名称=="APTT"),"项目结果"]=31#<2019-06-24
  # 2198476
  
  wuhu3[[84]][which(wuhu3[[84]]$项目名称=="APTT"),"项目结果"]=250#<2019-06-24

  
  wuhu3[[85]][which(wuhu3[[85]]$项目名称=="APTT"),"项目结果"]=250#<2019-06-24
  wuhu3[[85]][which(wuhu3[[85]]$项目名称=="PT(a)"),"项目结果"]=1.5#>2018-01-08
  wuhu3[[85]][which(wuhu3[[85]]$项目名称=="PT(r)"),"项目结果"]=40#>2018-01-08
  wuhu3[[85]][which(wuhu3[[85]]$项目名称=="PT(INR)"),"项目结果"]=39.99#>2018-01-08
  wuhu3[[85]][which(wuhu3[[85]]$项目名称=="PT(s)"),"项目结果"]=250#<2019-6-24
  wuhu3[[85]][which(wuhu3[[85]]$项目名称=="AT-III"),"项目结果"]=0#<2019-6-24
  
  wuhu3[[86]][which(wuhu3[[86]]$项目名称=="GA%"),"项目结果"]=13.8#median

  wuhu3[[87]][which(wuhu3[[87]]$项目名称=="HFLC%"),"项目结果"]=0.2
  wuhu3[[87]][which(wuhu3[[87]]$项目名称=="IG%"),"项目结果"]=3.0
  wuhu3[[87]][which(wuhu3[[87]]$项目名称=="P-LCR"),"项目结果"]=46.9
  wuhu3[[87]][which(wuhu3[[87]]$项目名称=="NRBC%"),"项目结果"]=9.2
  
  for(i in 88:93){
    wuhu3[[i]][which(wuhu3[[i]]$项目名称=="TT"),"项目结果"]=200#<2019-06-24
  }
  
  wuhu3[[94]][which(wuhu3[[94]]$项目名称=="APTT"),"项目结果"]=18.0#<2019-06-24
  # 2215550
  
  wuhu3[[95]][which(wuhu3[[95]]$项目名称=="APTT"),"项目结果"]=250#<2019-06-24

  wuhu3[[96]][which(wuhu3[[96]]$项目名称=="HFLC%"),"项目结果"]=0.15
  wuhu3[[96]][which(wuhu3[[96]]$项目名称=="IG%"),"项目结果"]=0.45
  wuhu3[[96]][which(wuhu3[[96]]$项目名称=="NRBC%"),"项目结果"]=0.15
  
  wuhu3[[97]][which(wuhu3[[97]]$项目名称=="HFLC%"),"项目结果"]=3.2
  wuhu3[[97]][which(wuhu3[[97]]$项目名称=="IG%"),"项目结果"]=3.9
  wuhu3[[97]][which(wuhu3[[97]]$项目名称=="P-LCR"),"项目结果"]=47.4
  wuhu3[[97]][which(wuhu3[[97]]$项目名称=="NRBC%"),"项目结果"]=1.5

  wuhu3[[98]][which(wuhu3[[98]]$项目名称=="APTT"),"项目结果"]=200#>=2019-06-24
  wuhu3[[98]][which(wuhu3[[98]]$项目名称=="TT"),"项目结果"]=250#>=2019-6-24
  wuhu3[[98]][which(wuhu3[[98]]$项目名称=="PT(a)"),"项目结果"]=1.5#>2018-01-08
  wuhu3[[98]][which(wuhu3[[98]]$项目名称=="PT(r)"),"项目结果"]=40#>2018-01-08
  wuhu3[[98]][which(wuhu3[[98]]$项目名称=="PT(INR)"),"项目结果"]=39.99#>2018-01-08
  wuhu3[[98]][which(wuhu3[[98]]$项目名称=="PT(s)"),"项目结果"]=450#>=2019-6-24
  wuhu3[[98]][which(wuhu3[[98]]$项目名称=="FIB"),"项目结果"]=0.15#>=2019-6-24

  wuhu3[[99]][which(wuhu3[[99]]$项目名称=="APTT"),"项目结果"]=250#<2019-06-24

  wuhu3[[100]][which(wuhu3[[100]]$项目名称=="TT"),"项目结果"]=12.1#>=2019-6-24

  wuhu3[[101]][which(wuhu3[[101]]$项目名称=="TT"),"项目结果"]=250#<2019-06-24

}
resultALL_no_na=rbindlist(c(list(resultALL_no_na),wuhu3))
names(resultALL_no_na)==names(resultALL)
resultALL=subset(resultALL_no_na,select = colnames(resultALL))
resultALL$项目结果=as.numeric(resultALL$项目结果)




djh=(unique(resultALL$就诊卡号))#13497
bah=(unique(resultALL$病案号))#13492
o2=dplyr::distinct(resultALL[,c("就诊卡号","病案号")])
o2=subset(o2,!is.na(o2$病案号))
table(is.na(unique(o2$病案号)))
table(is.na(o2$病案号))
resultALL=merge(resultALL,o2,by="就诊卡号",all.x = T)
table(is.na(resultALL$病案号.y))
table(is.na(resultALL$病案号.x))
resultALL$病案号.x=resultALL$病案号.y
resultALL$病案号.y=NULL
colnames(resultALL)[colnames(resultALL)=="病案号.x"]="病案号"

resultALL$报告时间 <- lubridate::ymd_hms(resultALL$报告时间,tz=Sys.timezone(location = TRUE))
resultALL$标本采集时间 <- lubridate::ymd_hms(resultALL$标本采集时间,tz=Sys.timezone(location = TRUE))

resultALL %>%
  as_tibble() %>%
  mutate (across(where(is.character),toupper)) ->resultALL
resultALL=as.data.frame(resultALL)
resultALL$参考值=gsub(" - | -|- ","-",resultALL$参考值)

save(resultALL,file="resultALL.Rdata")
# load(file="resultALL.Rdata")
# 11128295
resultALL=distinct(resultALL)#少了35个
# 11128259
# distinct 删除

table(resultALL$项目名称,useNA="ifany")
table(is.na(resultALL$病案号),useNA="ifany")
table(is.na(resultALL$就诊卡号),useNA="ifany")

if(T){
  #注意先数据变换倍数变换 再改单位 参比 再处理项目中文注释 再处理项目名称 
  # 11128294
  resultALL=subset(resultALL,resultALL$项目名称!="镜检")
  resultALL=subset(resultALL,resultALL$项目名称!="SF-LDH")
  resultALL=subset(resultALL,resultALL$项目名称!="SFT-LDH")
  resultALL=subset(resultALL,resultALL$项目名称!="SAMPLENO(SPE)")
  
  # 11128283
  resultALL[which(resultALL$项目名称=="ABO血型"),"项目名称"]="ABO"
  resultALL[which(resultALL$项目名称=="ABO"),"项目中文注释"]="ABO血型(正反定型)"
  resultALL=subset(resultALL,项目名称!="ABO")
  # 11126091
  resultALL[which(resultALL$项目名称=="RH血型"),"项目名称"]="RH"
  resultALL[which(resultALL$项目名称=="RH"),"项目中文注释"]="RH血型"
  resultALL=subset(resultALL,项目名称!="RH")
  #n11123899
 resultALL[which(resultALL$项目名称=="MCH"),"项目中文注释"]="平均红细胞血红蛋白量"
  
 resultALL[which(resultALL$项目名称=="MCHC"),"项目中文注释"]="平均红细胞血红蛋白浓度"
  
 resultALL[which(resultALL$项目名称=="MCV"),"项目中文注释"]="平均红细胞容积"
  
 resultALL[which(resultALL$项目名称=="GGT"&resultALL$项目中文注释=="R-谷氨酰基转移酶"),"项目中文注释"]="R-谷氨酰基转肽酶"
  
 resultALL[which(resultALL$项目名称=="G-ALB"&resultALL$项目中文注释=="白蛋白"),"项目中文注释"]="糖化试剂盒测的白蛋白"
 
 resultALL[which(resultALL$项目名称=="WBC"),"单位"]="×10^9/L"
  
 resultALL[which(resultALL$项目名称=="MONO%"),"项目中文注释"]="单核细胞百分数"
  
 resultALL[which(resultALL$项目名称=="MONO#"),"单位"]="×10^9/L"
  
 resultALL[which(resultALL$项目名称=="TG"),"项目中文注释"]="甘油三酯"
  
 resultALL[which(resultALL$项目名称=="HFLC%"),"单位"]="%"
 resultALL[which(resultALL$项目名称=="HFLC%"),"项目中文注释"]="高荧光细胞百分比"
  
 resultALL[which(resultALL$项目名称=="RBC"),"单位"]="×10^12/L"
  
 resultALL[which(resultALL$项目名称=="CK-MB"&resultALL$项目中文注释=="肌酸激酶同工酶"&resultALL$参考值=="7月25日"),"参考值"]="7-25"
  
 resultALL[which(resultALL$项目名称=="LYMPH#"),"单位"]="×10^9/L"
  
 resultALL[which(resultALL$项目名称=="LYMPH%"&resultALL$项目中文注释=="淋巴细胞比率"),"项目中文注释"]="淋巴细胞百分数"
  
 resultALL[which(resultALL$项目名称=="CSF-MTP"&resultALL$单位=="MG/L"),"项目结果"]=
   resultALL[which(resultALL$项目名称=="CSF-MTP"&resultALL$单位=="MG/L"),"项目结果"]/10
 resultALL[which(resultALL$项目名称=="CSF-MTP"),"单位"]="MG/DL"
 resultALL[which(resultALL$项目名称=="CSF-MTP"&resultALL$单位=="MG/DL"),"参考值"]="15.0-45.0"
  
  #连着转
 resultALL[which(resultALL$项目名称=="PT(R)"),"单位"]="RATIO"
  #
 resultALL[which(resultALL$项目名称=="PT(INR)"&resultALL$项目中文注释=="凝血酶原时间国际标准化比值"),"单位"]="INR"
  
 resultALL[which(resultALL$项目名称=="BASO%"&resultALL$项目中文注释=="嗜碱性粒细胞比率"),"项目中文注释"]="嗜碱性粒细胞百分数"
  
 resultALL[which(resultALL$项目名称=="BASO#"),"单位"]="×10^9/L"
  
 resultALL[which(resultALL$项目名称=="LIP"),"项目名称"]="LPS"
  
 resultALL[which(resultALL$项目名称=="EO%"),"项目名称"]="EOS%"
 resultALL[which(resultALL$项目名称=="EOS%"&resultALL$项目中文注释=="嗜酸性粒细胞比率"),"项目中文注释"]="嗜酸性粒细胞百分数"
  
 resultALL[which(resultALL$项目名称=="EO#"),"项目名称"]="EOS#"
 resultALL[which(resultALL$项目名称=="EOS#"&resultALL$项目中文注释=="嗜酸性粒细胞绝对值"),"单位"]="×10^9/L"
  
 resultALL[which(resultALL$项目名称=="GA"&resultALL$项目中文注释=="糖化白蛋白"&resultALL$单位=="G/DL"),"项目中文注释"]="糖化白蛋白_定量"
 resultALL[which(resultALL$项目名称=="GA"&resultALL$项目中文注释=="糖化白蛋白_定量"&resultALL$单位=="G/DL"),"项目名称"]="GA_DL"
 resultALL[which(resultALL$项目名称=="GA%"&resultALL$项目中文注释=="糖化白蛋白百分比"),"项目中文注释"]="糖化白蛋白"
 resultALL[which(resultALL$项目名称=="GA%"&resultALL$项目中文注释=="糖化白蛋白百分比"&resultALL$单位=="%"),"项目中文注释"]="糖化白蛋白"
 resultALL[which(resultALL$项目名称=="GA%"&resultALL$项目中文注释=="糖化白蛋白"),"项目名称"]="GA"  
 resultALL[which(resultALL$项目名称=="GA"&resultALL$项目中文注释=="糖化白蛋白"&resultALL$单位=="%"),"参考值"]="11.0-17.0"  
  
 resultALL[which(resultALL$项目名称=="RET‰"&resultALL$项目中文注释=="网织红细胞千分比"&resultALL$单位=="‰"),"项目结果"]=
   resultALL[which(resultALL$项目名称=="RET‰"&resultALL$项目中文注释=="网织红细胞千分比"&resultALL$单位=="‰"),"项目结果"]/10
 resultALL[which(resultALL$项目名称=="RET‰"&resultALL$项目中文注释=="网织红细胞千分比"&resultALL$单位=="‰"&resultALL$参考值=="6.5-16.9"),"参考值"]="0.65-1.69"
 resultALL[which(resultALL$项目名称=="RET‰"&resultALL$项目中文注释=="网织红细胞千分比"&resultALL$单位=="‰"&resultALL$参考值=="6.4-15.2"),"参考值"]="0.64-1.52"
 resultALL[which(resultALL$项目名称=="RET‰"&resultALL$项目中文注释=="网织红细胞千分比"&resultALL$单位=="‰"),"项目中文注释"]="网织红细胞百分比"
 resultALL[which(resultALL$项目名称=="RET‰"&resultALL$单位=="‰"),"单位"]="%"
 resultALL[which(resultALL$项目名称=="RET‰"),"项目名称"]="RET%"
 resultALL[which(resultALL$项目名称=="RET%"&resultALL$项目中文注释=="网织红细胞百分比"&resultALL$单位=="%"&resultALL$患者性别=="男"&is.na(resultALL$参考值)),"参考值"]="0.65-1.69"
 resultALL[which(resultALL$项目名称=="RET%"&resultALL$项目中文注释=="网织红细胞百分比"&resultALL$单位=="%"&resultALL$患者性别=="女"&is.na(resultALL$参考值)),"参考值"]="0.64-1.52"  
  
 resultALL[which(resultALL$项目名称=="HGB"),"项目名称"]="HB"
  
 resultALL[which(resultALL$项目名称=="PLT"&resultALL$单位=="G/L"),"单位"]="×10^9/L"
 resultALL[which(resultALL$项目名称=="PLT"&resultALL$参考值=="100.0-300.0"),"参考值"]="100-300"
  
 resultALL[which(resultALL$项目名称=="NRBC%"&resultALL$项目中文注释=="有核红细胞百分比"),"单位"]="%"
 resultALL[which(resultALL$项目名称=="NRBC%"&resultALL$项目中文注释=="有核红细胞百分比"),"项目中文注释"]="有核红细胞百分数"
  
 resultALL[which(resultALL$项目名称=="IG%"&resultALL$项目中文注释=="幼稚粒细胞百分比"),"单位"]="%"
 resultALL[which(resultALL$项目名称=="IG%"&resultALL$项目中文注释=="幼稚粒细胞百分比"),"项目中文注释"]="幼稚粒细胞百分数"
  
 resultALL[which(resultALL$项目名称=="CRE"&resultALL$参考值=="45.0-84.0"),"参考值"]="45-84"
  
 resultALL[which(resultALL$项目名称=="LPA"&resultALL$项目中文注释=="脂蛋白(A)"&resultALL$单位=="NMOL/L"),"项目结果"]= resultALL[which(resultALL$项目名称=="LPA"&resultALL$项目中文注释=="脂蛋白(A)"&resultALL$单位=="NMOL/L"),"项目结果"]/2.5
  
 resultALL[which(resultALL$项目名称=="LPA"&resultALL$项目中文注释=="脂蛋白(A)"),"单位"]="MG/DL"
 resultALL[which(resultALL$项目名称=="LPA"&resultALL$项目中文注释=="脂蛋白(A)"),"参考值"]="<30.0"
  
 resultALL[which(resultALL$项目名称=="NEUT%"&resultALL$项目中文注释=="中性粒细胞比率"),"单位"]="%"
 resultALL[which(resultALL$项目名称=="NEUT%"&resultALL$项目中文注释=="中性粒细胞比率"),"项目中文注释"]="中性粒细胞百分数"
  
 resultALL[which(resultALL$项目名称=="NEUT#"&resultALL$项目中文注释=="中性粒细胞绝对值"),"单位"]="×10^9/L"
  
 resultALL[which(resultALL$项目名称=="有核红细胞"&resultALL$单位=="个/100WBC"),"项目结果"]=
   resultALL[which(resultALL$项目名称=="有核红细胞"&resultALL$单位=="个/100WBC"),"项目结果"]/100
 resultALL[which(resultALL$项目名称=="有核红细胞"),"项目名称"]="NRBC"
 resultALL[which(resultALL$项目名称=="NRBC"),"项目中文注释"]="有核红细胞"
 resultALL[which(resultALL$项目名称=="NRBC"),"单位"]="×10^9/L"  
 resultALL[which(resultALL$项目名称=="AFU"&resultALL$项目中文注释=="血清A-L-岩藻苷酶"),"项目中文注释"]="血清A-L-岩藻糖苷酶"
  
 resultALL[which(resultALL$项目名称=="BUN"),"项目名称"]="UREA"
 resultALL[which(resultALL$项目名称=="UREA"),"项目中文注释"]="尿素"

 resultALL=distinct(resultALL)
 save(resultALL,file="resultALL02.Rdata")
 # load(file="resultALL02.Rdata")

 
  #参考值处理   1000088391  生化整体错了
  o1= resultALL[which(resultALL$就诊卡号==1000088391&resultALL$标本采集时间=="2020-11-20 23:14:00"),]
  o2= resultALL[which(resultALL$就诊卡号==1000088391&resultALL$标本采集时间=="2020-11-09 07:47:31"),]
  resultALL=subset(resultALL,!(resultALL$就诊卡号==1000088391&resultALL$标本采集时间=="2020-11-20 23:14:00"))
  for(x in o2$项目名称){
o1[which(o1$就诊卡号==1000088391&o1$标本采集时间=="2020-11-20 23:14:00"&o1$项目名称==x),"参考值"]=
o2[which(o2$就诊卡号==1000088391&o2$标本采集时间=="2020-11-09 07:47:31"&o2$项目名称==x),"参考值"]
  }
  resultALL=rbind(resultALL,o1)
  #
  resultALL[which(resultALL$就诊卡号==1703846&resultALL$标本采集时间=="2017-01-09 09:12:57"&resultALL$项目名称=="ALP"),"参考值"]="45.0-135.0"
  #outlier 处理 1000088391
 resultALL[which(resultALL$就诊卡号==1000088391&resultALL$标本采集时间=="2020-11-20 23:14:00"&resultALL$项目名称=="ALP"),"参考值"]="45.0-135.0"
  #
 resultALL[which(resultALL$就诊卡号==1000106600&resultALL$标本采集时间=="2020-08-18 06:07:26"&resultALL$项目名称=="PALB"),"参考值"]="20.0-43.0"
  # 2056850
 resultALL[which(resultALL$就诊卡号==2056850&resultALL$标本采集时间=="2018-01-08 06:40:16"&resultALL$项目名称=="PT(A)"),"参考值"]="70-130"
 resultALL[which(resultALL$就诊卡号==2056850&resultALL$标本采集时间=="2018-01-08 06:40:16"&resultALL$项目名称=="PT(INR)"),"参考值"]="0.85-1.15"
  # UREA 没有问题的不需要调整 
  # CRE没有问题
  # IBIL
 resultALL[which(resultALL$就诊卡号==1000057152&resultALL$标本采集时间>"2020-08-07"&resultALL$标本采集时间<"2020-08-08"&resultALL$项目名称=="IBIL"),"参考值"]="<= 17.0"
 resultALL[which(resultALL$就诊卡号==82695&resultALL$标本采集时间>"2020-08-07"&resultALL$标本采集时间<"2020-08-08"&resultALL$项目名称=="IBIL"),"参考值"]="<= 22.0"

  #
 resultALL[which(resultALL$就诊卡号%in%c(1000057152,82695,1000106757,1000184401)&resultALL$标本采集时间>="2020-08-07"&resultALL$标本采集时间<"2021-01-04"&resultALL$项目名称=="GLU"&resultALL$参考值=="3.89-6.38"),"参考值"]="3.9-6.1"
  # TBIL
 resultALL[which(resultALL$标本采集时间>="2018-01-08"&resultALL$标本采集时间<"2018-01-09"&resultALL$项目名称=="TBIL"),"参考值"]="2-21"
 resultALL[which(resultALL$就诊卡号%in%c(82695,1000184401)&resultALL$标本采集时间>"2020-08-07"&resultALL$标本采集时间<"2021-01-04"&resultALL$项目名称=="TBIL"),"参考值"] ="<= 26.0"
 resultALL[which(resultALL$就诊卡号%in%c(1000057152,1000106757,1000050017)&resultALL$标本采集时间>"2020-08-07"&resultALL$标本采集时间<"2021-01-04"&resultALL$项目名称=="TBIL"),"参考值"] ="<= 21.0"
  #
 resultALL[which(resultALL$就诊卡号%in%c(1000057152,82695,1000106757,1000050017,1000184401)&resultALL$标本采集时间>"2020-08-07"&resultALL$标本采集时间<"2021-01-04"&resultALL$项目名称=="DBIL"&resultALL$参考值=="0-5.1"),"参考值"]="<= 4.0"
 #

 resultALL[is.na(resultALL$报告时间),]$报告时间 = resultALL[is.na(resultALL$报告时间),]$标本采集时间+60*60*3
 resultALL = arrange(resultALL, resultALL$标本采集时间, resultALL$报告时间)
 resultALL_semester = subset(resultALL, select = c(就诊卡号, 标本采集时间, 项目名称))
 resultALL_semester = duplicated(resultALL_semester)
 # manual check 后注释
 resultALL[resultALL_semester, ]$标本采集时间 = resultALL[resultALL_semester, ]$标本采集时间+1
 
 resultALL[which(resultALL$项目名称=="GLU"&resultALL$标本采集时间=="2021-06-03 07:27:03"),"项目结果"]=5.81#只弄了一次
 resultALL[which(resultALL$项目名称=="GGT"&resultALL$标本采集时间=="2021-05-12 09:10:53"),"项目结果"]=744.2
 
 #科室原因造成的
 #日间治疗	张春英
 #泌外门诊	张春英
 #妇科门诊	梁颖
 #腹内门诊	梁颖   
 
 #CRE没有问题
 # resultALL=distinct(resultALL)
 # resultALL=subset(resultALL,!is.na(resultALL$项目结果))
 resultALL_remain=subset(resultALL,!resultALL$项目名称 %in% c("P-LCR","TBA","RDW-CV","RDW-SD","PDW","MPV","PCT","PLT","ALB","G","A/G","A/G/ALB(SPE)","ALB(SPE)","Γ-G(SPE)","Β2-G(SPE)","Β1-G(SPE)","Α2-G(SPE)","Α1-G(SPE)"))
 # PDW RDW-CV RDW-SD

 
 resultALL_xxb=subset(resultALL,resultALL$项目名称 %in% c("P-LCR","TBA","RDW-CV","RDW-SD","PDW","MPV","PCT","PLT"))
 
 for(i in unique(resultALL_xxb[which(resultALL_xxb$项目名称=="P-LCR"&is.na(resultALL_xxb$项目结果)),]$就诊卡号)){
   resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="P-LCR"&is.na(resultALL_xxb$项目结果)),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="P-LCR"),]$项目结果,na.rm = T)}
 resultALL_xxb[which(resultALL_xxb$项目名称=="P-LCR"&is.na(resultALL_xxb$项目结果)),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$项目名称=="P-LCR"),]$项目结果,na.rm = T)
 
 
 for(i in unique(resultALL_xxb[which(resultALL_xxb$项目名称=="TBA"&is.na(resultALL_xxb$项目结果)),]$就诊卡号)){resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="TBA"&is.na(resultALL_xxb$项目结果)),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="TBA"),]$项目结果,na.rm = T)}

 for(i in unique(resultALL_xxb[which(resultALL_xxb$项目名称=="RDW-CV"&is.na(resultALL_xxb$项目结果)),]$就诊卡号)){resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="RDW-CV"&is.na(resultALL_xxb$项目结果)),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="RDW-CV"),]$项目结果,na.rm = T)}
 resultALL_xxb[which(resultALL_xxb$项目名称=="RDW-CV"&is.na(resultALL_xxb$项目结果)&resultALL_xxb$标本采集时间<"2020-11-01"),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$项目名称=="RDW-CV"&resultALL_xxb$标本采集时间<"2020-11-01"),]$项目结果,na.rm = T)
 resultALL_xxb[which(resultALL_xxb$项目名称=="RDW-CV"&is.na(resultALL_xxb$项目结果)&resultALL_xxb$标本采集时间>="2020-11-01"),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$项目名称=="RDW-CV"&resultALL_xxb$标本采集时间>="2020-11-01"),]$项目结果,na.rm = T)
#
 for(i in unique(resultALL_xxb[which(resultALL_xxb$项目名称=="RDW-SD"&is.na(resultALL_xxb$项目结果)),]$就诊卡号)){resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="RDW-SD"&is.na(resultALL_xxb$项目结果)),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="RDW-SD"),]$项目结果,na.rm = T)}
 resultALL_xxb[which(resultALL_xxb$项目名称=="RDW-SD"&is.na(resultALL_xxb$项目结果)&resultALL_xxb$标本采集时间<"2020-11-01"),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$项目名称=="RDW-SD"&resultALL_xxb$标本采集时间<"2020-11-01"),]$项目结果,na.rm = T)
 
 resultALL_xxb[which(resultALL_xxb$项目名称=="RDW-SD"&is.na(resultALL_xxb$项目结果)&resultALL_xxb$标本采集时间>="2020-11-01"),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$项目名称=="RDW-SD"&resultALL_xxb$标本采集时间>="2020-11-01"),]$项目结果,na.rm = T)
 #
 for(i in unique(resultALL_xxb[which(resultALL_xxb$项目名称=="PDW"&is.na(resultALL_xxb$项目结果)),]$就诊卡号)){resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="PDW"&is.na(resultALL_xxb$项目结果)),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="PDW"),]$项目结果,na.rm = T)}
 resultALL_xxb[which(resultALL_xxb$项目名称=="PDW"&is.na(resultALL_xxb$项目结果)&resultALL_xxb$标本采集时间<"2020-11-01"),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$项目名称=="PDW"&resultALL_xxb$标本采集时间<"2020-11-01"),]$项目结果,na.rm = T)
 resultALL_xxb[which(resultALL_xxb$项目名称=="PDW"&is.na(resultALL_xxb$项目结果)&resultALL_xxb$标本采集时间>="2020-11-01"),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$项目名称=="PDW"&resultALL_xxb$标本采集时间>="2020-11-01"),]$项目结果,na.rm = T)
 #
 for(i in unique(resultALL_xxb[which(resultALL_xxb$项目名称=="MPV"&is.na(resultALL_xxb$项目结果)),]$就诊卡号)){resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="MPV"&is.na(resultALL_xxb$项目结果)),]$项目结果=median(resultALL_xxb[which(resultALL_xxb$就诊卡号==i&resultALL_xxb$项目名称=="MPV"),]$项目结果,na.rm = T)}
 resultALL_xxb[which(resultALL_xxb$项目名称=="MPV"&is.na(resultALL_xxb$项目结果)),]$项目结果=10.3
 
 for(i in unique(resultALL_xxb[which(resultALL_xxb$项目名称=="PCT"&is.na(resultALL_xxb$项目结果)),]$ID_time)){resultALL_xxb[which(resultALL_xxb$ID_time==i&resultALL_xxb$项目名称=="PCT"&is.na(resultALL_xxb$项目结果)),]$项目结果=resultALL_xxb[which(resultALL_xxb$ID_time==i&resultALL_xxb$项目名称=="MPV"),]$项目结果*resultALL_xxb[which(resultALL_xxb$ID_time==i&resultALL_xxb$项目名称=="PLT"),]$项目结果/1000000}
 resultALL_remain=rbind(resultALL_remain,resultALL_xxb)
 # 10532394
 resultALL_calculate_before=subset(resultALL,项目名称 %in% c("ALB","G","A/G","A/G/ALB(SPE)","ALB(SPE)","Γ-G(SPE)","Β2-G(SPE)","Β1-G(SPE)","Α2-G(SPE)","Α1-G(SPE)"))
 #591505
 resultALL_calculate_before=unite(resultALL_calculate_before,col="ID_time2","就诊卡号","标本采集时间",sep = "_",remove = F)
 resultALL_calculate_before=split(resultALL_calculate_before,resultALL_calculate_before$ID_time2)
 table(unlist(lapply(resultALL_calculate_before,function(x)nrow(x)>10)))
 cl <- makeCluster(4)
 registerDoParallel(cl)
 resultALL_calculate_before=foreach(x=1:length(resultALL_calculate_before))%dopar%{
   x=dplyr::arrange(resultALL_calculate_before[[x]],报告时间)
   x$报告时间=x$报告时间[1]
   x=subset(x,select=colnames(resultALL_remain))
   return(x)
 }
 stopCluster(cl)
 resultALL_calculate_before=rbindlist(resultALL_calculate_before)
 # # 有三种情况
 # table((is.na(resultALL_calculate["PCT",])),(is.na(resultALL_calculate["MPV",])),(is.na(resultALL_calculate["PLT",])))
 # # 188188 全不缺
 # table((!is.na(resultALL_calculate["PCT",]))&(!is.na(resultALL_calculate["MPV",]))&(!is.na(resultALL_calculate["PLT",])))
 # # 137683 全缺
 # table((is.na(resultALL_calculate["PCT",]))&(is.na(resultALL_calculate["MPV",]))&(is.na(resultALL_calculate["PLT",])))
 # # 2422 就一个plt算不了哈哈哈
 # table((is.na(resultALL_calculate["PCT",]))&(is.na(resultALL_calculate["MPV",]))&(!is.na(resultALL_calculate["PLT",])))
 
 resultALL_calculate=resultALL_calculate_before
 # resultALL_calculate=distinct(resultALL_calculate[,c("项目名称","就诊卡号", "标本采集时间")])
 resultALL_calculate = tidyr::pivot_wider(
   resultALL_calculate,
   id_cols = "项目名称",
   names_from = c("就诊卡号", "标本采集时间"),
   names_prefix = "",
   names_sep = "_",
   names_glue = NULL,
   names_sort = FALSE,
   names_repair = "check_unique",
   values_from = "项目结果",
   values_fill = NA,
   values_fn = NULL
 )
 table(apply(resultALL_calculate,1,function(x) length(unlist(x))))
 row_names = resultALL_calculate[[1]]
 resultALL_calculate[[1]] = NULL
 col_names = colnames(resultALL_calculate)
 resultALL_calculate = as.matrix(resultALL_calculate)
 rownames(resultALL_calculate) = row_names
 colnames(resultALL_calculate) = col_names
 resultALL_calculate=as.data.frame(resultALL_calculate)
 resultALL_calculate["ALB","1700540_2015-08-24 07:22:14"] <- 43.1
 resultALL_calculate["ALB","1574621_2015-08-28 07:24:27"] <- 43.1
 resultALL_calculate["ALB","1708730_2015-08-30 07:53:10"] <- 37.0
 resultALL_calculate["ALB","1708730_2015-08-31 07:07:07"] <- 37.2
 resultALL_calculate["ALB","1708625_2015-08-30 07:53:12"] <- 36.6
 resultALL_calculate["ALB","1708625_2015-08-31 07:07:09"] <- 37.7
 resultALL_calculate["ALB","1441805_2015-08-31 07:07:08"] <- 39.6
 resultALL_calculate["ALB","1793799_2017-03-03 08:29:35"] <- 26.1
 #
 resultALL_calculate["TP","1378068_2016-10-20 07:27:39"] <- 74.9
 resultALL_calculate["G","1378068_2016-10-20 07:27:39"] <- 39.5
 resultALL_calculate["TP","1378068_2016-12-27 07:43:25"] <- 78.7
 resultALL_calculate["G","1378068_2016-12-27 07:43:25"] <- 42.7
 resultALL_calculate["TP","1378068_2016-09-22 07:34:40"] <- 75.5
 resultALL_calculate["G","1378068_2016-09-22 07:34:40"] <- 37.9
 
 resultALL_calculate["G","2205685_2019-03-05 06:23:12"] <- 21.9
 #
 resultALL_calculate["G","1865920_2021-09-16 07:02:06"] <- 31.7
#
 resultALL_calculate["A/G",]= resultALL_calculate["ALB",]/resultALL_calculate["G",]
 resultALL_calculate["A/G/ALB(SPE)",]= resultALL_calculate["A/G",]*(100-resultALL_calculate["ALB(SPE)",])/resultALL_calculate["ALB(SPE)",]
 
 # r3=lapply(resultALL_calculate,function(x) (!is.na(x[1]))&is.na(x[2]) )
 # r3=lapply(resultALL_calculate,function(x) (!is.na(x[3]))&is.na(x[2]) )
 # r3=lapply(resultALL_calculate,function(x) (!is.na(x[4]))&is.na(x[2]) )
 # r3=lapply(resultALL_calculate,function(x) (!is.na(x[5]))&is.na(x[2]) )
 # yyy=subset(resultALL_calculate,select=unlist(r3))
 # dput(names(yyy))
 # c("1342540_2014-07-21 06:37:00", "1213587_2014-07-21 06:38:00", 
 #   "1303428_2014-07-21 06:38:00", "1329817_2014-07-21 06:38:00", 
 #   "2133970_2018-05-21 11:21:43")
 # #前面是真没有啊
 # 
 # table(resultALL[which(resultALL$项目名称=="G-ALB"),"项目中文注释"])
 
 
 resultALL_calculate=dplyr::mutate(as.data.frame(resultALL_calculate),"项目名称"=rownames(resultALL_calculate))
 # resultALL_calculate=tidyr::as_tibble(resultALL_calculate)
 resultALL_calculate=tidyr::pivot_longer(
   data=resultALL_calculate,
   cols=!项目名称,#everthing
   names_to = c("就诊卡号","标本采集时间"),
   names_prefix = NULL,
   names_sep = "_",
   names_pattern = NULL,
   # names_ptypes = list(),
   # names_transform = list(),
   names_repair = "check_unique",
   values_to = "项目结果",
   values_drop_na = T,
   # values_ptypes = list(),
   # values_transform = list()
 )
 resultALL_calculate$标本采集时间=lubridate::ymd_hms((resultALL_calculate$标本采集时间),tz=Sys.timezone(location = TRUE))
 resultALL_calculate$就诊卡号=as.numeric(resultALL_calculate$就诊卡号)
 nrow(resultALL_calculate)#1042158
 # resultALL_calculate=merge(resultALL_calculate,resultALL_calculate_before,by=c("项目名称","就诊卡号","标本采集时间"),all.x=T)
 # resultALL_calculate$项目结果.y=NULL
 # colnames(resultALL_calculate)[colnames(resultALL_calculate)=="项目结果.x"]="项目结果"
 # resultALL_calculate=subset(resultALL_calculate,select=colnames(resultALL_calculate_before))
 # table(is.na(resultALL_calculate$项目中文注释))
 # nrow(resultALL_calculate)
 #TG
 # resultALL_calculate$项目名称=="TP"&
 # resultALL_calculate_median=resultALL_calculate[which(is.na(resultALL_calculate$项目中文注释)),c("就诊卡号","标本采集时间","项目名称","项目结果")]
 resultALL_calculate_median=merge(resultALL_calculate,resultALL_calculate_before,by=c("就诊卡号","标本采集时间"),all.x = T)
 resultALL_calculate_median$项目名称.y=NULL
 resultALL_calculate_median$项目结果.y=NULL
 colnames(resultALL_calculate_median)[colnames(resultALL_calculate_median)=="项目结果.x"]="项目结果"
 colnames(resultALL_calculate_median)[colnames(resultALL_calculate_median)=="项目名称.x"]="项目名称"

 # resultALL_calculate_median=resultALL_calculate_median2

for(i in 1:3){
 resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="A/G"&resultALL_calculate_median$标本采集时间>="2016-11-19"),c("项目中文注释", "单位", "参考值")[i]]=c("白蛋白/球蛋白比值",NA,"1.2-2.4")[i]
 }
 
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="A/G"&resultALL_calculate_median$标本采集时间<"2016-11-19"),c("项目中文注释", "单位", "参考值")[i]]=c("白蛋白/球蛋白比值",NA,"1.0-2.5")[i]
 }
 
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="ALB"&resultALL_calculate_median$标本采集时间>="2016-11-19"),c("项目中文注释", "单位", "参考值")[i]]=c("白蛋白","G/L","40.0-55.0")[i]
 }
 
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="ALB"&resultALL_calculate_median$标本采集时间<"2016-11-19"),c("项目中文注释", "单位", "参考值")[i]]=c("白蛋白","G/L","35-51")[i]
 }

 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="A/G/ALB(SPE)"&resultALL_calculate_median$标本采集时间>="2016-11-19"),c("项目中文注释", "单位", "参考值")[i]]=c("白蛋白/球蛋白比值与蛋白电泳白蛋白比值",NA,"batch2")[i]
 }
 
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="A/G/ALB(SPE)"&resultALL_calculate_median$标本采集时间<"2016-11-19"),c("项目中文注释", "单位", "参考值")[i]]=c("白蛋白/球蛋白比值与蛋白电泳白蛋白比值",NA,"batch1")[i]
 } 
 
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="TP"&resultALL_calculate_median$标本采集时间>="2016-11-19"),c("项目中文注释", "单位", "参考值")[i]]=c("总蛋白","G/L","65.0-85.0")[i]
 }
 
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="TP"&resultALL_calculate_median$标本采集时间<"2016-11-19"),c("项目中文注释", "单位", "参考值")[i]]=c("总蛋白","G/L","60-85")[i]
 } 

for(i in 1:3){
  resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="G"),c("项目中文注释", "单位", "参考值")[i]]=c("球蛋白","G/L","20-40")[i]
}


 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="ALB(SPE)"),c("项目中文注释", "单位", "参考值")[i]]=c("白蛋白(毛细管电泳法)","%","52.0-62.8")[i]
 }
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="Γ-G(SPE)"),c("项目中文注释", "单位", "参考值")[i]]=c("Γ球蛋白(毛细管电泳法)","%","13.1-23.3")[i]
 }
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="Α1-G(SPE)"),c("项目中文注释", "单位", "参考值")[i]]=c("Α1球蛋白(毛细管电泳法)","%","3.1-4.6")[i]
 }
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="Α2-G(SPE)"),c("项目中文注释", "单位", "参考值")[i]]=c("Α2球蛋白(毛细管电泳法)","%","7.0-11.1")[i]
 }
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="Β1-G(SPE)"),c("项目中文注释", "单位", "参考值")[i]]=c("Β1球蛋白(毛细管电泳法)","%","5.3-7.8")[i]
 }
 for(i in 1:3){
   resultALL_calculate_median[which(resultALL_calculate_median$项目名称=="Β2-G(SPE)"),c("项目中文注释", "单位", "参考值")[i]]=c("Β1球蛋白(毛细管电泳法)","%","3.3-6.4")[i]
 }

table(resultALL_calculate_median$项目名称)
 #标本采集时间>"2020-08-07"
 # &is.na(resultALL_calculate$项目中文注释)
 resultALL_calculate_median=distinct(resultALL_calculate_median)
 table(as.data.frame(table(resultALL_calculate_median$ID_time,deparse.level = 1))$Freq)
  wuhuhuhuhuhuhuhu=as.data.frame(table(resultALL_calculate_median$ID_time,deparse.level = 1))
wuhuhuhuhuhuhuhu=subset(wuhuhuhuhuhuhuhu,Freq>10)
resultALL_calculate_median=subset(resultALL_calculate_median,select = colnames(resultALL_remain))
 
 # resultALL_calculate[which(resultALL_calculate$项目名称=="A/G/ALB(SPE)"&is.na(resultALL_calculate$项目中文注释)&resultALL_calculate$标本采集时间>=ymd("2016-11-19",tz=Sys.timezone(location = TRUE))),]
 # 
 # resultALL_calculate[which(resultALL_calculate$项目名称=="A/G/ALB(SPE)"&is.na(resultALL_calculate$项目中文注释)&resultALL_calculate$标本采集时间<ymd("2016-11-19",tz=Sys.timezone(location = TRUE))),]
 
 resultALL=rbind(resultALL_remain,resultALL_calculate_median)
 
 resultALL=subset(resultALL,select =c("申请科室", "患者姓名", "患者性别", "患者年龄", "就诊卡号","病案号", "诊断", "标本采集时间", "专业名称", "报告时间", "项目名称","项目中文注释", "项目结果", "单位", "参考值"))


 # 0.6 1.4
 # 计算插补的max
 # if(F){
 #   # table(is.na(resultALL$项目结果))  
 #   # Hmisc::describe()
 #   # pastecs::
 #   #   psych::describe()
 #   # psych::describe(,na.rm = T)
 #   
 #   head(sort(resultALL[resultALL$项目名称=="AT-III"&resultALL$标本采集时间>="2019-06-24",]$项目结果,decreasing = F),n=1L,na.rm = T)[1]/1.5
 #   head(sort(resultALL[resultALL$项目名称=="AT-III"&resultALL$标本采集时间<"2019-06-24",]$项目结果,decreasing = F),n=1L,na.rm = T)[1]/1.5
 #   
 #   head(sort(resultALL[resultALL$项目名称=="TT"&resultALL$标本采集时间>="2019-06-24",]$项目结果,decreasing = T),n=1L,na.rm = T)[1]*1.5
 #   head(sort(resultALL[resultALL$项目名称=="TT"&resultALL$标本采集时间<"2019-06-24",]$项目结果,decreasing = T),n=1L,na.rm = T)[1]*1.5
 #   
 #   head(sort(resultALL[resultALL$项目名称=="D-D"&resultALL$标本采集时间<"2019-06-24",]$项目结果,decreasing = T),n=20L,na.rm = T)[1]*1.5
 #   
 #   head(sort(resultALL[resultALL$项目名称=="PT(S)"&resultALL$标本采集时间>="2019-06-24",]$项目结果,decreasing = T),n=1L,na.rm = T)[1]*1.5
 #   head(sort(resultALL[resultALL$项目名称=="PT(S)"&resultALL$标本采集时间<"2019-06-24",]$项目结果,decreasing = T),n=1L,na.rm = T)[1]*1.5
 #   
 #   head(sort(resultALL[resultALL$项目名称=="PT(INR)"&resultALL$标本采集时间>="2018-01-08",]$项目结果,decreasing = T),n=1L,na.rm = T)[1]*1.5
 #   head(sort(resultALL[resultALL$项目名称=="PT(INR)"&resultALL$标本采集时间<"2018-01-08",]$项目结果,decreasing = T),n=1L,na.rm = T)[1]*1.5
 #   
 #   head(sort(resultALL[resultALL$项目名称=="PT(A)"&resultALL$标本采集时间>="2018-01-08",]$项目结果,decreasing = F),n=1L,na.rm = T)[1]/1.5
 #   head(sort(resultALL[resultALL$项目名称=="PT(A)"&resultALL$标本采集时间<"2018-01-08",]$项目结果,decreasing = F),n=1L,na.rm = T)[1]/1.5
 #   
 #   head(sort(resultALL[resultALL$项目名称=="FIB"&resultALL$标本采集时间>="2019-06-24",]$项目结果,decreasing = F),n=1L,na.rm = T)[1]/1.5
 #   head(sort(resultALL[resultALL$项目名称=="FIB"&resultALL$标本采集时间<"2019-06-24",]$项目结果,decreasing = F),n=1L,na.rm = T)[1]/1.5
 #   
 #   head(sort(resultALL[resultALL$项目名称=="APTT"&resultALL$标本采集时间>="2019-06-24",]$项目结果,decreasing = T),n=1L,na.rm = T)[1]*1.5
 #   head(sort(resultALL[resultALL$项目名称=="APTT"&resultALL$标本采集时间<"2019-06-24",]$项目结果,decreasing = T),n=1L,na.rm = T)[1]*1.5
 #   
 #   head(sort(resultALL[resultALL$项目名称=="FDP",]$项目结果,decreasing = T),n=1L,na.rm = T)[1]*1.5#1500
 #   
 #   head(sort(resultALL[resultALL$项目名称=="LPA",]$项目结果,decreasing = T),n=1L,na.rm = T)[1]*1.5
 #   
 #   median(resultALL[resultALL$项目名称=="APTT"&resultALL$标本采集时间<"2019-06-24",]$项目结果,na.rm = T)
 #   median(resultALL[resultALL$项目名称=="TT"&resultALL$标本采集时间<"2019-06-24",]$项目结果,na.rm = T)
 #   
 #   
 #   median(resultALL[resultALL$项目名称=="D-D"&resultALL$标本采集时间<"2019-06-24",]$项目结果,na.rm = T)
 #   
 #   median(resultALL[resultALL$项目名称=="FIB"&resultALL$标本采集时间<"2019-06-24",]$项目结果,na.rm = T)   
 #   median(resultALL[resultALL$项目名称=="AT-III"&resultALL$标本采集时间<"2019-06-24",]$项目结果,na.rm = T)
 #   
 #   
 #   head(sort(resultALL[resultALL$项目名称=="PT(INR)",]$项目结果,decreasing = T),n=10L,na.rm = T)
 #   head(sort(resultALL[resultALL$项目名称=="NRBC%",]$项目结果,decreasing = T),n=10L,na.rm = T)
 #   median(resultALL[resultALL$项目名称=="NRBC%",]$项目结果,na.rm=T)
 #   
 #   
 #   median(resultALL[resultALL$项目名称=="FIB"&resultALL$标本采集时间<"2019-06-24",]$项目结果,na.rm = T)
 #   
 #   median(resultALL[resultALL$项目名称=="PDW"&resultALL$就诊卡号==2247369,]$项目结果,na.rm = T)
 #   median(resultALL[resultALL$项目名称=="P-LCR"&resultALL$就诊卡号==2247369,]$项目结果,na.rm = T)
 #   median(resultALL[resultALL$项目名称=="IG%"&resultALL$就诊卡号==2247369,]$项目结果,na.rm = T)
 #   median(resultALL[resultALL$项目名称=="NRBC%"&resultALL$就诊卡号==2247369,]$项目结果,na.rm = T)
 #   
 #   median(resultALL[resultALL$项目名称=="GA%",]$项目结果,na.rm = T)#
 #   
 # }
 
 #下面处理错误的参考值，在上方必须将处理为唯一
 save(resultALL,file="resultALL after.Rdata")
 # load(file="resultALL after.Rdata")
}



resultALLarrange=dplyr::arrange(resultALL, resultALL$标本采集时间)
resultALLarrange=split(resultALLarrange,resultALLarrange$就诊卡号)
resultALLarrange[[1]]$标本采集时间[nrow(resultALLarrange[[1]])]
resultALLarrangeSEQ=lapply(resultALLarrange,function(x) difftime(x$标本采集时间[nrow(x)], x$标本采集时间[1],units = "weeks"))
table(resultALLarrangeSEQ>=10)
grep("NULL",names(resultALLarrangeSEQ))
table(is.na(names(resultALLarrangeSEQ)))
oooooooo=names(resultALLarrangeSEQ)[resultALLarrangeSEQ>=10]
resultALL3month=resultALL[resultALL$就诊卡号%in%oooooooo,]
oooooooo=dplyr::distinct(resultALL3month[,c("就诊卡号","病案号")])
table(is.na(oooooooo$病案号))
save(resultALL,resultALL3month,oooooooo,resultALLarrange,file="jyk14210831.Rdata")
# load(file="jyk14210831.Rdata")
resultALL3month_yx=subset(resultALL3month,grepl("胰",resultALL3month$诊断))

resultALL3month_yx=subset(resultALL3month_yx,resultALL3month_yx$标本采集时间<=as.Date("2019-12-31"))
resultALL3month_yx=subset(resultALL3month_yx,resultALL3month_yx$申请科室%in%c("胰胃外科1","胰胃外科2","腹外科3"))
unique(resultALL3month_yx$就诊卡号)
openxlsx::write.xlsx(resultALL3month_yx,file = "resultALL3month_yx.xlsx")

resultALL3month_yx_2021=subset(resultALL3month_yx,resultALL3month_yx$标本采集时间>as.Date("2019-12-31"))
resultALL3month_yx_2021=subset(resultALL3month_yx_2021,resultALL3month_yx_2021$申请科室%in%c("胰胃外科1","胰胃外科2","腹外科3"))
unique(resultALL3month_yx_2021$就诊卡号)
openxlsx::write.xlsx(resultALL3month_yx_2021,file = "resultALL3month_yx_20200101.xlsx")

gc1417_dsk=readxl::read_excel("C:/Users/zxh/Desktop/R/近全远端胃筛选/1417远端匹配检验科.xlsx",sheet=1)
gc1417_dsk=gc1417_dsk[gc1417_dsk$登记号 %in%oooooooo$就诊卡号,]
oooooooo=as.data.frame(oooooooo)
writexl::write_xlsx(gc1417_dsk,path = "C:/Users/zxh/Desktop/R/近全远端胃筛选/匹配上的有三个月以上的检验科数据.xlsx")
writexl::write_xlsx(as.data.frame(oooooooo),path = "C:/Users/zxh/Desktop/R/近全远端胃筛选/匹配上的有三个月以上的检验数据登记号和病案号.xlsx")

resultALL_enoughnum = resultALL
resultALL_enoughnum$标本采集时间 = quarter(resultALL_enoughnum$标本采集时间, with_year = TRUE)
if(F){
pdf('Figure A.Box before del.pdf',width = 20,height = 400)
box = resultALL_enoughnum
box$项目结果 = log2(resultALL_enoughnum$项目结果+1)
if (require('ggpubr')) {
  library(ggpubr)
  p <-ggboxplot(box,x = "标本采集时间",y = "项目结果",#group="项目名称",# color = "项目名称",
      font.label = list(size = 1, color = "black"),
      outlier.shape = 1#,palette=c(1,2,3),order=c("No changes","Shallow Deletion","Deep Deletion")#"#00AFBB","#E7B800",
    ) + labs(title = "", x = " ", y = "Expression")#+#tag="ENSG00000000003"
  p + facet_grid(rows = vars(box$项目名称), margins = TRUE) + stat_compare_means(method = "kruskal.test", label.y =2)
}
dev.off()
}
resultALL_enoughnum = as.data.frame(table(resultALL_enoughnum$标本采集时间, resultALL_enoughnum$项目名称))
resultALL_enoughnum = pivot_wider(
  resultALL_enoughnum,
  id_cols = "Var2",
  names_from = "Var1",
  names_prefix = "",
  names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "check_unique",
  values_from = "Freq",
  # values_fill = NULL,
  values_fn = NULL
)

resultALL_enoughnum[resultALL_enoughnum == 0] <- NA
resultALL_enoughnum = as.data.frame(resultALL_enoughnum)
rownames(resultALL_enoughnum) = resultALL_enoughnum$Var2
resultALL_enoughnum$Var2 = NULL
resultALL_enoughnum = resultALL_enoughnum[rowSums(resultALL_enoughnum, na.rm = T) > 10000, ]
sort(apply(resultALL_enoughnum, 1, function(x)  sum(is.na(x)) / ncol(resultALL_enoughnum)))
resultALL_enoughnum = resultALL_enoughnum[apply(resultALL_enoughnum, 1, function(x) sum(is.na(x)) / ncol(resultALL_enoughnum)) < 1 / 3, ]
resultALL_enoughnum %>% filter_at(colnames(.), any_vars(. < 500)) -> resultALL_enoughnum_needhebing

resultALL = subset(resultALL, resultALL$项目名称 %in% rownames(resultALL_enoughnum))
resultALL = distinct(resultALL)
if(F){
pdf('Figure A.Boxpl.pdf', width = 20, height = 400)
box = resultALL
box$标本采集时间 = quarter(box$标本采集时间, with_year = TRUE)
box$项目结果 = log2(resultALL$项目结果+1)
if (require('ggpubr')) {
  library(ggpubr)
  p <-
    ggboxplot(
      box,
      x = "标本采集时间",
      y = "项目结果",
      #group="项目名称",# color = "项目名称",
      font.label = list(size = 1, color = "black"),
      outlier.shape = 1#,palette=c(1,2,3),order=c("No changes","Shallow Deletion","Deep Deletion")#"#00AFBB","#E7B800",
    ) + labs(title = "", x = " ", y = "Expression")#+#tag="ENSG00000000003"
  p + facet_grid(rows = vars(box$项目名称), margins = TRUE) + stat_compare_means(method = "kruskal.test", label.y =2)
}
dev.off()
}
table(resultALL$项目名称)

resultALL_szsj=resultALL
resultALL_szsj %>% mutate (across(where(is.POSIXct),  ~ as.POSIXlt(., tz =Sys.timezone(location = TRUE)))) -> resultALL_szsj
resultALL_szsj %>% mutate (across(where(is.POSIXlt),  ~ as_date(., tz = Sys.timezone(location = TRUE)))) -> resultALL_szsj
resultALL_szsj = subset(resultALL_szsj, select = c("就诊卡号", "患者姓名", "病案号", "标本采集时间"))
resultALL_szsj = dplyr::distinct(resultALL_szsj)
table(is.na(resultALL_szsj$标本采集时间))
# ascend
resultALL_szsj_accend = dplyr::arrange(resultALL_szsj, resultALL_szsj$标本采集时间)
kanquexiang = resultALL_szsj_accend[resultALL_szsj_accend$就诊卡号  %in% resultALL_szsj_accend[is.na(resultALL_szsj_accend$病案号),]$就诊卡号,]#所有没有病案的都是这个看的
resultALL_szsj_accend = resultALL_szsj_accend[!duplicated(resultALL_szsj_accend$就诊卡号),]

table(duplicated(resultALL_szsj_accend$就诊卡号))
table(is.na(resultALL_szsj_accend$病案号))#407
# resultALL_szsj_accend = resultALL_szsj_accend[!is.na(resultALL_szsj_accend$病案号), ]#不要没住院
resultALL_szsj_accend_formerge = resultALL_szsj_accend[, c("就诊卡号","病案号", "标本采集时间")]
colnames(resultALL_szsj_accend_formerge) = c("djh","bah", "szsj")
######
# resultALL[resultALL_semester,]$报告时间=resultALL[resultALL_semester,]$报告时间+1

table(apply(resultALL_enoughnum, 1, function(x) sum(is.na(x))))
# 0  1  4  9
# 24 55  5  2
resultALL_semester = resultALL[, c("就诊卡号", "标本采集时间", "报告时间", "项目名称", "项目结果")]

resultALL_list <- list()
resultALL_list[[1]] = subset(resultALL_semester,resultALL_semester$项目名称 %in% rownames(resultALL_enoughnum[apply(resultALL_enoughnum, 1, function(x)sum(is.na(x))) == 0, ]))#没有问题直接插补
resultALL_list[[2]] = subset(resultALL_semester,resultALL_semester$项目名称 %in% rownames(resultALL_enoughnum[apply(resultALL_enoughnum, 1, function(x)sum(is.na(x))) == 1, ]))#要处理
resultALL_list[[3]] = subset(resultALL_semester,resultALL_semester$项目名称 %in% rownames(resultALL_enoughnum[apply(resultALL_enoughnum, 1, function(x)sum(is.na(x))) == 4, ]))#要处理
resultALL_list[[4]] = subset(resultALL_semester,resultALL_semester$项目名称 %in% rownames(resultALL_enoughnum[apply(resultALL_enoughnum, 1, function(x)sum(is.na(x))) == 9, ]))#没有问题直接插补
resultALL_list = lapply(1:4, function(x) {
  y = tidyr::pivot_wider(
    resultALL_list[[x]] ,
    id_cols = "项目名称",
    names_from = c("就诊卡号", "标本采集时间"),
    names_prefix = "",
    names_sep = "_",
    names_glue = NULL,
    names_sort = FALSE,
    names_repair = "check_unique",
    values_from = "项目结果",
    values_fill = NULL,
    values_fn = NULL
  )
  row_names = y[[1]]
  y[[1]] = NULL
  col_names = colnames(y)
  y = as.matrix(y)
  rownames(y) = row_names
  colnames(y) = col_names
  y
})
# check dup reason
if (F) {
  wula <- list()
  for (i in 1:4) {
    wula[[i]] = apply(resultALL_list[[i]], 1, function(x)
      length(unlist(x)) > length(x))
    print(table(wula[[i]]))
  }
}
lapply(1:4, function(x) table(is.na(resultALL_list[[x]])))

resultALL_list[[5]] = subset(
  resultALL_semester,
  项目名称  %in% c(
    "PT(R)","PT(A)","APTT","PT(INR)","FIB","PT(S)","TT","AT-III","FDP","D-D"
  )
)
need_process = c(rownames(resultALL_list[[2]]), rownames(resultALL_list[[3]]))
need_process = need_process[!need_process %in% c("PT(R)", "PT(A)", "APTT", "PT(INR)", "FIB", "PT(S)", "TT", "AT-III", "FDP", "D-D")]
resultALL_list[[6]] = subset(resultALL_semester, 项目名称  %in% need_process)
resultALL_list[[2]] = {
  y = tidyr::pivot_wider(
    resultALL_list[[5]] ,
    id_cols = "项目名称",
    names_from = c("就诊卡号", "标本采集时间"),
    names_prefix = "",
    names_sep = "_",
    names_glue = NULL,
    names_sort = FALSE,
    names_repair = "check_unique",
    values_from = "项目结果",
    values_fill = NULL,
    values_fn = NULL
  )
  row_names = y[[1]]
  y[[1]] = NULL
  col_names = colnames(y)
  y = as.matrix(y)
  rownames(y) = row_names
  colnames(y) = col_names
  y
}
resultALL_list[[3]] = {
  y = tidyr::pivot_wider(
    resultALL_list[[6]] ,
    id_cols = "项目名称",
    names_from = c("就诊卡号", "标本采集时间"),
    names_prefix = "",
    names_sep = "_",
    names_glue = NULL,
    names_sort = FALSE,
    names_repair = "check_unique",
    values_from = "项目结果",
    values_fill = NULL,
    values_fn = NULL
  )
  row_names = y[[1]]
  y[[1]] = NULL
  col_names = colnames(y)
  y = as.matrix(y)
  rownames(y) = row_names
  colnames(y) = col_names
  y
}
lapply(1:length(resultALL_list),function(x) table(is.na(resultALL_list[[x]])))

# resultALL_list[[5]]=NULL
resultALL_list[[5]] = resultALL_list[[3]][apply(resultALL_list[[3]], 1, function(x)
  sum(is.na(x))) <= 5000, ]
resultALL_list[[5]] = resultALL_list[[5]][, apply(resultALL_list[[5]], 2, function(x)
  sum(is.na(x)) != length(x))]
resultALL_list[[6]] = resultALL_list[[3]][apply(resultALL_list[[3]], 1, function(x)
  sum(is.na(x))) > 5000, ]
resultALL_list[[6]] = resultALL_list[[6]][, apply(resultALL_list[[6]], 2, function(x)
  sum(is.na(x)) != length(x))]
resultALL_list[[3]] = NULL
lapply(1:length(resultALL_list), function(x)  table(is.na(resultALL_list[[x]])))
  

  resultALL_sexdiff=resultALL[, c("项目名称","患者性别",#"标本采集时间",
                                  "项目中文注释", "单位", "参考值")]
  resultALL_sexdiff %>%group_by_all() %>% mutate(count = n()) -> resultALL_sexdiff#非常好不懂顺序
  resultALL_sexdiff$标本采集时间=resultALL$标本采集时间
  # resultALL_sexdiff=subset(resultALL_sexdiff,resultALL_sexdiff$项目名称%in%resultALL_adjust$项目名称)
  resultALL_sexdiff=arrange(resultALL_sexdiff,desc(resultALL_sexdiff$标本采集时间))
  resultALL_sexdiff$标本采集时间=as_date(resultALL_sexdiff$标本采集时间,tz=Sys.timezone(location = TRUE))
  resultALL_sexdiff=distinct(resultALL_sexdiff)
  #质控看有没把日期分错的
  if(F){
    pdf('Figure A.Boxplot of 重叠1.pdf',width=400,height=8)
    p=ggplot(resultALL_sexdiff, aes(fill=参考值,x=标本采集时间)) + geom_bar(position="fill")+facet_grid(cols = vars(resultALL_sexdiff$项目名称), margins = TRUE)+theme(legend.position="none")#,group=参考值
    p
    dev.off()

  }
  resultALL_sexdiff1 = resultALL_sexdiff[, c("项目名称","标本采集时间","参考值")]
  resultALL_sexdiff1 = resultALL_sexdiff1[!duplicated(resultALL_sexdiff1[, c("项目名称","参考值")],fromLast = T), ]
  resultALL_sexdiff1=distinct(resultALL_sexdiff1)
  colnames(resultALL_sexdiff1)[2]= "最早标本采集时间"

  resultALL_sexdiff = resultALL_sexdiff[!duplicated(resultALL_sexdiff[, c("项目名称","患者性别", "项目中文注释", "单位", "参考值")]), ]
  resultALL_sexdiff = arrange(resultALL_sexdiff, resultALL_sexdiff$患者性别)
  resultALL_sexdiff = resultALL_sexdiff[!duplicated(resultALL_sexdiff[, c("项目名称", "项目中文注释", "单位", "参考值")]), ]
  resultALL_sexdiff=arrange(resultALL_sexdiff,resultALL_sexdiff$项目名称)
  resultALL_sexdiff=distinct(resultALL_sexdiff)
  resultALL_sexdiff=merge(resultALL_sexdiff,resultALL_sexdiff1,by=c("项目名称","参考值"))

  resultALL_sexdiff1=split(resultALL_sexdiff,resultALL_sexdiff$项目名称)
  
  for(i in 1:length(resultALL_sexdiff1)){
    resultALL_sexdiff1[[i]]%>% group_by(标本采集时间)%>% mutate(batch1=group_indices()) ->     resultALL_sexdiff1[[i]]
  }
  
  if(T){
  #ALP
  # 其实2020-7-15 也碰到这个问题
  # #   <=2016-11-18 15岁我看不出来
  # <=2016-11-18（包）<=15 任何 0-750
  # <=2016-11-18（包）>=16 任何 40-150
  # #   <=2016-11-19 <= <=2017-01-09
  # #男
  # >=2016-11-19  age <=12 男 35.0-300.0 (重复)
  # 2016-11-19<=date<=2017-1-9  13<=age<=17  男   0.0-500.0
  # 2017-01-10<=date<=2021-09-22  13<=age<=17 男 40.0-390.0
  # 2016-11-19<=date<=2017-1-9 (包) age>=18 男 45.0-135.0
  # 2017-01-10<=date<=2021-09-22 (包) age>=18 男 45.0-125.0
  # 
  # #女
  # >=2016-11-19  age <=12 女 35.0-300.0 (重复)
  # >=2016-11-19 (包) 13<=age<=17 女 35.0-187.0
  # >=2016-11-19 (包) 18<=age<=49 女 35.0-100.0
  # >=2016-11-19 (包) 50<=age 女 50.0-135.0
  resultALL_sexdiff1[["ALP"]]$batch1=1
  resultALL_sexdiff1[["ALP"]][resultALL_sexdiff1[["ALP"]]$参考值%in%c("40-150","0-750"),"batch1"]=2
  # CRE
  # <=2016-11-18（包）任何 男 59-104
  # <=2016-11-18（包）任何 女 45-84
  # 老人组一致
  # >=2016-11-19 (包) >=60 男 57.0-111.0
  # >=2016-11-19 (包) >=60 女 41.0-81.0
  # 幼年组一致 在 2020-7-15 到 2020-8-6 的某天似乎被暂停
  # >=2016-11-19 (包) age<=19 男 59-104
  # >=2016-11-19 (包) age<=19 女 45-84
  # >2020-7-15 (包) age<=19 男 57.0-97.0
  # >2020-7-15 (包) age<=19 女 45-84
  # 青年组
  # >=2016-11-19 (包) 18<=age<=59 男 57.0-97.0
  # >=2016-11-19 (包) 18<=age<=59 女 41.0-73.0
  #按照男女分开，年龄分列1-3-2  
  resultALL_sexdiff1[["CRE"]]$batch1=2
  # resultALL_sexdiff1[["CRE"]][resultALL_sexdiff1[["CRE"]]$参考值%in%c("40-150","0-750"),"batch1"]=2
  # resultALL_sexdiff1[["CRE"]][resultALL_sexdiff1[["CRE"]]$参考值=="41.0-73.0",]
    
  # resultALL_sexdiff1[["CRE"]][nrow(resultALL_sexdiff1[["CRE"]])+1,]=  c(
  CRE1=resultALL_sexdiff1[["CRE"]][resultALL_sexdiff1[["CRE"]]$参考值=="59-104",]
  CRE1["最早标本采集时间"]=as_date("2016-11-19",tz=Sys.timezone(location = TRUE))
  CRE2=resultALL_sexdiff1[["CRE"]][resultALL_sexdiff1[["CRE"]]$参考值=="45-84",]
  CRE2["最早标本采集时间"]=as_date("2016-11-19",tz=Sys.timezone(location = TRUE))
  resultALL_sexdiff1[["CRE"]][resultALL_sexdiff1[["CRE"]]$参考值%in%c("59-104","45-84"),"标本采集时间"]=as_date("2016-11-18",tz=Sys.timezone(location = TRUE))
  resultALL_sexdiff1[["CRE"]][resultALL_sexdiff1[["CRE"]]$参考值%in%c("59-104","45-84"),"batch1"]=1
  resultALL_sexdiff1[["CRE"]]=rbind(  resultALL_sexdiff1[["CRE"]],CRE1,CRE2)
  # UREA
  # <=2016-11-18（包）任何 任何  2.14-7.14 
  # 老人组一致
  # >=2016-11-19 (包) >=60 男  3.6-9.5 
  # >=2016-11-19 (包) >=60 女  3.1-8.8 
  # 幼年组一致 在 2020-7-15 到 2020-8-6 的某天似乎被暂停
  # >=2016-11-19 (包) age<=19 任何  2.14-7.14 
  # >2020-7-15 (包) age<=19 男  3.1-8.0
  # >2020-7-15 (包) age<=19 女  2.6-7.5
  # 青年组
  # >=2016-11-19 (包) 18<=age<=59 男  3.1-8.0 
  # >=2016-11-19 (包) 18<=age<=59 女  2.6-7.5

  #HCY
  #0-60岁  0-15.0
  #>=61岁  0-20.0
  

  }
  resultALL_semester_batch=lapply(resultALL_sexdiff1,function(x)as.data.frame(x))
  resultALL_semester_batch=data.table::rbindlist(resultALL_semester_batch)
  length(unique(resultALL[is.na(resultALL$病案号),]$就诊卡号))
  openxlsx::write.xlsx(resultALL_semester_batch,file ="resultALL_adjust_ref_date.xlsx")#writexl 保存原本的信息
  
  table(sort(c(resultALL_semester_batch$标本采集时间,resultALL_semester_batch$最早标本采集时间)))
  
  save(resultALL,resultALL_sexdiff1,resultALL_list,resultALL_semester_batch,file="resultALL after1.Rdata")
  # load(file="resultALL after1.Rdata")

  resultALL_sexdiff2=lapply(resultALL_sexdiff1,function(x) max(x$batch1))
  resultALL_sexdiff2=tidyr::as_tibble(resultALL_sexdiff2)
  

  resultALL_semester_batch3_0=dput(colnames(resultALL_sexdiff2)[resultALL_sexdiff2[1,]==3])
  # 2018-01-07 2018-01-08 2020-08-06 2020-08-07 
  # c("TBIL")#直接在上面处理了
  # ALP TBIL UREA CRE

  
  # 
  #
  #
  dput(colnames(resultALL_sexdiff2)[resultALL_sexdiff2[1,]==2])
  # 2016-11-18 2016-11-19
  # 2018-01-07 2018-01-08
  # 2019-06-23 2019-06-24
  # 2020-08-06 2020-08-07
  # 2020-10-31 2020-11-01
  resultALL_semester_batch2=subset(resultALL_semester_batch,resultALL_semester_batch$项目名称 %in% c("A/G","A/G/ALB(SPE)", "ALB", "ALP", "ALT", "APOA1", "APOB", "APTT", "AST", "AT-III", "BASO#", "CA", "CHOL", "CK", "CK-MB", "CL", "CRE","D-D", "DBIL", "EOS#", "EOS%", "FE", "FIB", "GGT", "GLU", "HB","HBDH", "HCT", "HDL-CHO", "IBIL", "IGA", "IGG", "IGM", "K", "LDH","LDL-CHO", "LYMPH#", "LYMPH%", "MCHC", "MCV", "MG", "MONO#","MONO%", "NA", "NEUT#", "NEUT%", "PALB", "PCT", "PDW", "PHOS",    "PLT", "PT(A)", "PT(INR)", "PT(S)", "RBC", "RDW-CV", "RDW-SD","TCO2", "TG", "TP", "TRANSFE", "TT", "UREA", "URIC", "WBC", "Β2-MG"))
  resultALL_semester_batch2_1=unique(resultALL_semester_batch2[which(resultALL_semester_batch2$最早标本采集时间==as_date("2016-11-19")),"项目名称"])
  resultALL_semester_batch2_2=unique(resultALL_semester_batch2[which(resultALL_semester_batch2$最早标本采集时间==as_date("2018-01-08")),"项目名称"])
  resultALL_semester_batch2_3=unique(resultALL_semester_batch2[which(resultALL_semester_batch2$最早标本采集时间==as_date("2019-06-24")),"项目名称"])
  resultALL_semester_batch2_4=unique(resultALL_semester_batch2[which(resultALL_semester_batch2$最早标本采集时间==as_date("2020-08-07")),"项目名称"])
  resultALL_semester_batch2_5=unique(resultALL_semester_batch2[which(resultALL_semester_batch2$最早标本采集时间==as_date("2020-11-01")),"项目名称"])
  
  table(resultALL_semester_batch2$项目名称 %in%c(resultALL_semester_batch2_1$项目名称,resultALL_semester_batch2_2$项目名称,resultALL_semester_batch2_3$项目名称,resultALL_semester_batch2_4$项目名称,resultALL_semester_batch2_5$项目名称))
  #
  resultALL_semester_batch1_0=dput(colnames(resultALL_sexdiff2)[resultALL_sexdiff2[1,]==1])
  # resultALL_sexdiff1_1=subset(wula,wula$项目名称 %in% c("ADA", "ALB(SPE)", "BASO%", "CRP", "FDP", "G", "HCY", "LPA", "MCH", "MPV", "P-LCR", "PT(R)", "SOD", "TBA", "Α1-G(SPE)", "Α2-G(SPE)", "Β1-G(SPE)", "Β2-G(SPE)", "Γ-G(SPE)"))
  resultALL_semester_batch_ref=list(resultALL_semester_batch1_0,resultALL_semester_batch2_1$项目名称,resultALL_semester_batch2_2$项目名称,resultALL_semester_batch2_3$项目名称,resultALL_semester_batch2_4$项目名称,resultALL_semester_batch2_5$项目名称,resultALL_semester_batch3_0)

    # cl <- makeCluster(7)
    # registerDoParallel(cl)
    # foreach(i=1:7,.packages = c(library(dplyr),library(lubridate)))%dopar%{

# for(i in 1:7){
    # box=subset(resultALL,resultALL$项目名称=="TBIL")
    # 2016-11-18 2016-11-19
    # 2018-01-07 2018-01-08
    # 2019-06-23 2019-06-24
    # 2020-08-06 2020-08-07
    # 2020-10-31 2020-11-01
    # box=subset(resultALL,项目名称 %in% resultALL_semester_batch_ref[[i]])

if(F){
    pdf(file =paste0("Figure A. resultALL_semester_batch_ref.pdf"),height = 2000,width = 20)
    box=resultALL
    box=box %>% mutate(res = case_when(
      标本采集时间 %within% interval("1900-01-01","2016-11-18 23:59:59",tz=Sys.timezone(location = TRUE)) ~ "20140504 to 20161118",#1
      标本采集时间 %within% interval("2016-11-19","2018-01-07 23:59:59",tz=Sys.timezone(location = TRUE)) ~ "20161119 to 20180107",#2
      标本采集时间 %within% interval("2018-01-08","2019-06-23 23:59:59",tz=Sys.timezone(location = TRUE)) ~ "20180108 to 20190623",#3
      标本采集时间 %within% interval("2019-06-24","2020-08-06 23:59:59",tz=Sys.timezone(location = TRUE))  ~ "20190624 to 20200806",#4
      标本采集时间 %within% interval("2020-08-07","2020-10-31 23:59:59",tz=Sys.timezone(location = TRUE))  ~ "20200807 to 20201031",#5
      标本采集时间 %within% interval("2020-11-01","2500-01-01",tz=Sys.timezone(location = TRUE))  ~ "20201101 to 2021922",#6
      TRUE~"wrong"))
    box$项目结果=log2(box$项目结果+1)
    
    p <- ggpubr::ggboxplot(box, x = "res", y = "项目结果",color = "res",#group="项目名称", 
                   font.label = list(size = 1, color = "black"),outlier.shape = NA#,palette=c(1,2,3),order=c("No changes","Shallow Deletion","Deep Deletion")#"#00AFBB","#E7B800",
    )+labs(title="",x=" ", y = "log2 ( Laboratory results + 1)")+stat_compare_means(method = "kruskal.test", label.y =9)+stat_compare_means(aes(color="res"),label.y =seq(5, 9, length.out = 15),method = "wilcox.test", comparisons =combn(1:6, 2, FUN = list))+facet_grid(rows = vars(box$项目名称),space="free_y",margins = F)#+facet_wrap(~项目名称,nrow = 6)#+#tag="ENSG00000000003"
    p
    #主要问题在2020-08-06 2020-08-07，而前面比较相似
    dev.off()
}
    # deparse(substitute(box))
    # eval(parse(text = "box"))
    # deparse(substitute(box))

  
    #直接用rowmeans来插补
    # k <- which(is.na(wula), arr.ind=TRUE)
    # wula[k] <- rowMeans(wula, na.rm=TRUE)[k[,1]]

    # resultALL_list_NA_location <- list()
    # for(i in 1:length(resultALL_list)){
    #   resultALL_list_NA_location[[i]]  <- which(is.na(resultALL_list[[i]]), arr.ind=TRUE)
    #   resultALL_list[[i]][resultALL_list_NA_location[[i]]] <- rowMeans(resultALL_list[[i]], na.rm=TRUE)[resultALL_list_NA_location[[i]][,1]]
    #   # mod = model.matrix(~as.factor(cancer), data=pheno)
    # }
    # 
 
if(T){    
    #处理tbil as single
    resultALL_semester_batch_ref_TBIL=list(resultALL_semester_batch2_1$项目名,resultALL_semester_batch2_2$项目名称,resultALL_semester_batch2_3$项目名称,c(resultALL_semester_batch2_4$项目名称,resultALL_semester_batch3_0),resultALL_semester_batch2_5$项目名称)#resultALL_semester_batch1_0 这些不需要变化的
  wula <- list()
  wula=foreach(j=1:length(resultALL_semester_batch_ref_TBIL))%:%foreach(i=1:length(resultALL_list))%do%{subset(resultALL_list[[i]],rownames(resultALL_list[[i]]) %in% resultALL_semester_batch_ref_TBIL[[j]])}
  lapply(wula,lengths)
  # [[1]]
  # [1]       0       0       0 2090074  482150
  # 
  # [[2]]
  # [1]      0 113914      0      0      0
  # 
  # [[3]]
  # [1]      0 341742      0      0      0
  # 
  # [[4]]
  # [1]       0       0  124892  597164 1350020
  # 
  # [[5]]
  # [1] 3812200       0       0       0       0
  wula=unlist(wula,recursive = F)
  wula=wula[lengths(wula)>0]
  wula=lapply(wula,function(x) as.data.frame(x))
  lapply(wula,function(x) table(unlist(lapply(x,function(x) sum(is.na(x))/length(x)==1))))
  batch=lapply(wula,function(x) colnames(x)[unlist(lapply(x,function(x) sum(is.na(x))/length(x)!=1))])
  wula=lapply(1:length(wula),function(x) subset(wula[[x]],select=batch[[x]]))
  lapply(wula,function(x) table(unlist(lapply(x,function(x) sum(is.na(x))/length(x)==1))))
  #这个就效率很低不知道为啥
  #如果多对于1个那么需要改成 case when if else
  batch <- list()
  for (i in 1:2) {
    batch[[i]]=lapply(wula,function(x) as.numeric(ymd_hms(stringr::str_split(colnames(x),"_",simplify=T)[,2],tz=Sys.timezone(location = TRUE)) %within% interval("1900-01-01","2016-11-18 23:59:59",tz=Sys.timezone(location = TRUE))))[[i]]
  }
  for (i in 3) {
    batch[[i]]=lapply(wula,function(x) as.numeric(ymd_hms(stringr::str_split(colnames(x),"_",simplify=T)[,2],tz=Sys.timezone(location = TRUE)) %within% interval("1900-01-01","2018-01-07 23:59:59",tz=Sys.timezone(location = TRUE))))[[i]]
  }
  for (i in 4) {
    batch[[i]]=lapply(wula,function(x) as.numeric(ymd_hms(stringr::str_split(colnames(x),"_",simplify=T)[,2],tz=Sys.timezone(location = TRUE)) %within% interval("1900-01-01","2019-06-23 23:59:59",tz=Sys.timezone(location = TRUE))))[[i]]
  }
  for (i in 5:7) {
    batch[[i]]=lapply(wula,function(x) as.numeric(ymd_hms(stringr::str_split(colnames(x),"_",simplify=T)[,2],tz=Sys.timezone(location = TRUE)) %within% interval("1900-01-01","2020-08-06 23:59:59",tz=Sys.timezone(location = TRUE))))[[i]]
  }
  for (i in 8) {
    batch[[i]]=lapply(wula,function(x) as.numeric(ymd_hms(stringr::str_split(colnames(x),"_",simplify=T)[,2],tz=Sys.timezone(location = TRUE)) %within% interval("1900-01-01","2020-10-31 23:59:59",tz=Sys.timezone(location = TRUE))))[[i]]
  }
  # 先按照batch 日期来分
  # 再按照result_list来分类
  wula_NA_location <- list()
  cl <- makeCluster(8)
  registerDoParallel(cl)
  wula_new <- list()
  wula_new=foreach(i=1:8)%dopar%{
  wula_NA_location[[i]]  <- which(is.na(wula[[i]]), arr.ind=TRUE)
  if(length(wula_NA_location[[i]])!=0){
  wula[[i]][wula_NA_location[[i]]] <- rowMeans(wula[[i]], na.rm=TRUE)[wula_NA_location[[i]][,1]]
  }
  wula_new =sva::ComBat(dat =wula[[i]],batch = batch[[i]],mod=NULL,mean.only = TRUE, ref.batch=0)
  if(length(wula_NA_location[[i]])!=0){
  wula_new[wula_NA_location[[i]]] <- NA
  }
  wula_new=dplyr::mutate(as.data.frame(wula_new),"项目名称"=rownames(wula_new))
  wula_new=tidyr::as_tibble(wula_new)
  wula_new=tidyr::pivot_longer(
    data=wula_new,
    cols=!项目名称,#everthing
    names_to = c("就诊卡号","标本采集时间"),
    names_prefix = NULL,
    names_sep = "_",
    names_pattern = NULL,
    # names_ptypes = list(),
    # names_transform = list(),
    names_repair = "check_unique",
    values_to = "项目结果",
    values_drop_na = T,
    # values_ptypes = list(),
    # values_transform = list()
  )
  wula_new$标本采集时间=lubridate::ymd_hms((wula_new$标本采集时间),tz=Sys.timezone(location = TRUE))
  wula_new$就诊卡号=as.numeric(wula_new$就诊卡号)
  wula_new$项目结果=as.numeric(wula_new$项目结果)
  return(wula_new)
  }
  stopCluster(cl)
  wula_new=rbindlist(wula_new)
  wula_new=as.data.frame(wula_new)
  wula_no_needadjust=subset(resultALL,resultALL$项目名称%in%resultALL_semester_batch1_0,select=c("项目名称","就诊卡号","标本采集时间","项目结果"))
  nrow(resultALL)==(nrow(wula_no_needadjust)+nrow(wula_new))
  if(F){
    qwerty=subset(resultALL,!resultALL$项目名称%in%resultALL_semester_batch1_0,select=c("项目名称","就诊卡号","标本采集时间","项目结果"))
    qwerty=unite(qwerty,col="ID_time","项目名称","就诊卡号","标本采集时间",sep = "_",remove = F)
    wula_new123=unite(wula_new,col="ID_time","项目名称","就诊卡号","标本采集时间",sep = "_",remove = F)
    qwerty2=qwerty[!qwerty$ID_time %in%wula_new123$ID_time,]
    table(qwerty2$项目名称)
    }
  # GGT    GLU    PDW RDW-CV RDW-SD 
  # 1      1   2419    452    454 
  resultALL_ref_sva=rbind(wula_new,wula_no_needadjust)
  resultALL_ref_sva=merge(resultALL,resultALL_ref_sva,by = c("就诊卡号","标本采集时间","项目名称"))
  resultALL_ref_sva[which(resultALL_ref_sva$项目结果.y<=0),"项目结果.y"]=0
  resultALL_ref_sva_ori=resultALL_ref_sva
  resultALL_ref_sva$项目结果.x=NULL
  colnames(resultALL_ref_sva)[colnames(resultALL_ref_sva)=="项目结果.y"]="项目结果"
  # colnames(resultALL_ref_sva$项目结果.y)
  # $项目结果.y
  resultALL_ref_sva=subset(resultALL_ref_sva,select =c("申请科室", "患者姓名", "患者性别", "患者年龄", "就诊卡号","病案号", "诊断", "标本采集时间", "专业名称", "报告时间", "项目名称","项目中文注释", "项目结果", "单位", "参考值"))
  save(resultALL,resultALL_ref_sva_ori,resultALL_ref_sva,file="resultALL sva.Rdata")
  

  if(F){
  pdf(file =paste0("Figure A. resultALL sva.pdf"),height = 2000,width = 20)
  box=resultALL_ref_sva
  box=box %>% mutate(res = case_when(
    标本采集时间 %within% interval("1900-01-01","2016-11-18 23:59:59",tz=Sys.timezone(location = TRUE)) ~ "20140504 to 20161118",#1
    标本采集时间 %within% interval("2016-11-19","2018-01-07 23:59:59",tz=Sys.timezone(location = TRUE)) ~ "20161119 to 20180107",#2
    标本采集时间 %within% interval("2018-01-08","2019-06-23 23:59:59",tz=Sys.timezone(location = TRUE)) ~ "20180108 to 20190623",#3
    标本采集时间 %within% interval("2019-06-24","2020-08-06 23:59:59",tz=Sys.timezone(location = TRUE))  ~ "20190624 to 20200806",#4
    标本采集时间 %within% interval("2020-08-07","2020-10-31 23:59:59",tz=Sys.timezone(location = TRUE))  ~ "20200807 to 20201031",#5
    标本采集时间 %within% interval("2020-11-01","2500-01-01",tz=Sys.timezone(location = TRUE))  ~ "20201101 to 2021922",#6
    TRUE~"wrong"))
  box$项目结果.y=log2(box$项目结果.y+1)
  
  p <- ggpubr::ggboxplot(box, x = "res", y = "项目结果.y",color = "res",#group="项目名称", 
                         font.label = list(size = 1, color = "black"),outlier.shape = NA,order=c("20140504 to 20161118","20161119 to 20180107","20180108 to 20190623","20190624 to 20200806","20200807 to 20201031","20201101 to 2021922")#,palette=c(1,2,3)#"#00AFBB","#E7B800",
  )+labs(title="",x=" ", y = "log2 ( Laboratory results + 1)")+stat_compare_means(method = "kruskal.test", label.y =9)+stat_compare_means(aes(color="res"),label.y =seq(5, 9, length.out = 15),method = "wilcox.test", comparisons =combn(1:6, 2, FUN = list))+facet_grid(rows = vars(box$项目名称),space="free_y",margins = F)#+facet_wrap(~项目名称,nrow = 6)#+#tag="ENSG00000000003"
  p
  #主要问题在2020-08-06 2020-08-07，而前面比较相似
  dev.off()
  }
}  
    
########1 全胃

    JDW=read_excel(r"(C:\Users\zxh\Desktop\R\近端胃小篇论文作图表\1421全胃13号最终版(TSN修改2).xlsx)",sheet = 1,skip = 2)
    colnames(JDW)=stringr::str_split(colnames(JDW),"_",simplify = T)[,1]
    JDW=subset(JDW,wwww==1,select = !duplicated(colnames(JDW)))
    JDW=subset(JDW,select =!colnames(JDW)%in% c("CHOL","TG","LDL","HDL","Lpa","ApoB","ApoA1"))
    JDW$bah=as.numeric(JDW$bah)
    colnames(JDW)[duplicated(colnames(JDW))]
    JDW %>% mutate (across(where(is.POSIXct),~as.POSIXlt(., tz=Sys.timezone(location = TRUE)))) ->JDW
    JDW %>% mutate (across(where(is.POSIXlt),~as_date(., tz=Sys.timezone(location = TRUE)))) ->JDW
    resultALL_szsj_accend_formerge_havebah=resultALL_szsj_accend_formerge[!is.na(resultALL_szsj_accend_formerge$bah),]
    JDW=merge(JDW,resultALL_szsj_accend_formerge_havebah,by = "bah",all.x = T)
    # JDW[JDW$bah=="1144886","szsj"] <- "2014-06-09"#毛康京
    JDW$qzryts=JDW$rysj-JDW$szsj#确诊到入院总天数
    JDW$qzssts=JDW$ssrq-JDW$szsj#确诊到手术总天数
    if(F){
      # JDW=subset(JDW,wwww==1,select =c("bah", "name","rysj","cysj","ssrq","NAC"))
      writexl::write_xlsx(JDW,path = "JDW_maybewrong.xlsx")
      #拉出来需要重新录得
    }   
  
    
  JDW_replustime=read_excel(r"(C:\Users\zxh\Desktop\R\jianyanke\JDW_maybewrongdsk_20220212.xlsx)")
  colnames(JDW_replustime)=stringr::str_split(colnames(JDW_replustime),"_",simplify = T)[,1]
  dput(names(JDW_replustime))
  JDW_replustime=subset(JDW_replustime,select = c("bah","diffzd", "needcheck", "realfirstdate", "Lastneotransfusiondate","stopneooraldate", "specialdatenote1", "manualsqsjcalculate", "specialdatenote2"))
  JDW_replustime$realfirstdate=as.Date(JDW_replustime$realfirstdate,tz=Sys.timezone(location = TRUE))
  JDW_replustime$Lastneotransfusiondate=as.Date(JDW_replustime$Lastneotransfusiondate,tz=Sys.timezone(location = TRUE))
  JDW_replustime$stopneooraldate=as.Date(JDW_replustime$stopneooraldate,tz=Sys.timezone(location = TRUE))
  JDW_replustime$bah=as.numeric(JDW_replustime$bah)
  JDW_replustime=merge(JDW,JDW_replustime,by= "bah")
  table(JDW_replustime$needcheck)
  JDW_replustime[which(JDW_replustime$needcheck%in%c(1,2)),"szsj"] <-  JDW_replustime[which(JDW_replustime$needcheck%in%c(1,2)),"realfirstdate"]
  JDW_replustime$neoadjuvant=as.numeric(JDW_replustime$NAC==1|JDW_replustime$NAR==1)
  #可以手术的时间 dssj
  JDW_replustime$dssj=JDW_replustime$szsj
  JDW_replustime[which(JDW_replustime$NAC==1|JDW_replustime$NAR==1),"dssj"]=JDW_replustime[which(JDW_replustime$NAC==1|JDW_replustime$NAR==1),"stopneooraldate"]
  
  JDW_replustime$qzryts=JDW_replustime$rysj-JDW_replustime$szsj#确诊到入院总天数
  JDW_replustime$qzssts=JDW_replustime$ssrq-JDW_replustime$szsj#确诊到手术总天数
  table(  JDW_replustime$qzryts<0)
  table(  JDW_replustime$qzssts<0)

  # wwwww=JDW_replustime[JDW_replustime$dsts!=JDW_replustime$manualsqsjcalculate,]
  # write.csv(wwwww,"wwwww.csv")
  #院前天数（门诊收治到入院）
  JDW_replustime$yqts=JDW_replustime$rysj-JDW_replustime$dssj
  #术前时间
  JDW_replustime$sqts=JDW_replustime$ssrq-JDW_replustime$rysj
  #待术天数 = 院前天数 + 术前天数
  JDW_replustime$dsts=JDW_replustime$ssrq-JDW_replustime$dssj
  #术后天数（入院到手术）
  JDW_replustime$shts=JDW_replustime$cysj-JDW_replustime$ssrq
  #住院天数
  JDW_replustime$zyts=JDW_replustime$cysj-JDW_replustime$rysj

  JDW_replustime$After2020=ifelse(JDW_replustime$rysj>=as_date("2020-01-01"),1,0)
  # 下级医院到入院 直接手术
  JDW_replustime$xjyqts=JDW_replustime$yqts
  JDW_replustime[which(JDW_replustime$neoadjuvant==0),]$xjyqts=JDW_replustime[which(JDW_replustime$neoadjuvant==0),]$yqts+JDW_replustime[which(JDW_replustime$neoadjuvant==0),]$diffzd
  # 下级医院到手术 直接手术
  JDW_replustime$xjdsts=JDW_replustime$dsts
  JDW_replustime[which(JDW_replustime$neoadjuvant==0),]$xjdsts=JDW_replustime[which(JDW_replustime$neoadjuvant==0),]$dsts+JDW_replustime[which(JDW_replustime$neoadjuvant==0),]$diffzd
  # &JDW_replustime$neoadjuvant==0
  # JDW_replustime$diffzd[is.na(JDW_replustime$diffzd)] <- 0
  
  JDW_replustime <- JDW_replustime %>%
    mutate(across(where(is.difftime),as.numeric))
  quantile(JDW_replustime$sqts,probs = seq(0, 1, 0.1))
  JDW_replustime$sqts1 = cut(JDW_replustime$sqts,breaks = c(-Inf,4,Inf),right =T,ordered_result = T,labels = c("≤4",">4"))
  
  quantile(JDW_replustime$shts)
  JDW_replustime$shts1 = cut(JDW_replustime$shts,breaks = c(-Inf,14,Inf),right =T,ordered_result = T,labels = c("≤14",">14"))
  JDW_replustime$shts3=NA
  JDW_replustime$shts3=ifelse(JDW_replustime$shts>14,1,ifelse(JDW_replustime$shts<=14,0,"wrong"))
  JDW_replustime$shts3=as.numeric(JDW_replustime$shts3)
  
  quantile(JDW_replustime$zyts)
  JDW_replustime$zyts1 = cut(JDW_replustime$zyts,breaks = c(-Inf,14,Inf),right = T,ordered_result = T,labels = c("≤14",">14"))
  
  quantile(JDW_replustime$xjyqts)
  JDW_replustime$xjyqts1 = cut(JDW_replustime$xjyqts,breaks = c(-Inf,28,Inf),right =T,ordered_result = T,labels =c("≤28",">28"))
  
  JDW_replustime$xjyqts2 = cut(JDW_replustime$xjyqts,breaks = c(-Inf,30,Inf),right =T,ordered_result = T,labels =c("≤30",">30"))
  
  JDW_replustime$xjyqts3=NA
  JDW_replustime$xjyqts3=ifelse(JDW_replustime$xjyqts>30,1,ifelse(JDW_replustime$xjyqts<=30,0,"wrong"))
  JDW_replustime$xjyqts3=as.numeric(JDW_replustime$xjyqts3)
  
  quantile(JDW_replustime$qzryts)
  JDW_replustime$qzryts1 = cut(JDW_replustime$qzryts,breaks = c(-Inf,14,Inf),right =T,ordered_result = T,labels = c("≤14",">14"))
  quantile(JDW_replustime$qzssts)
  JDW_replustime$qzssts1 = cut(JDW_replustime$qzssts,breaks = c(-Inf,21,Inf),right =T,ordered_result = T,labels = c("≤21",">21"))
  
  quantile(JDW_replustime$xjdsts,probs = seq(0, 1, 0.1))
  quantile(JDW_replustime[which(JDW_replustime$neoadjuvant==1),]$xjdsts,probs = seq(0, 1, 0.1))
  quantile(JDW_replustime[which(JDW_replustime$neoadjuvant==0),]$xjdsts,probs = seq(0, 1, 0.1))
  
  JDW_replustime$xjdsts1 = cut(JDW_replustime$xjdsts,breaks = c(-Inf,28,Inf),right =T,ordered_result = T,labels =c("≤28",">28"))
  
  JDW_replustime$xjdsts2 = cut(JDW_replustime$xjdsts,breaks = c(-Inf,30,Inf),right =T,ordered_result = T,labels =c("≤30",">30"))
  
  JDW_replustime$xjdsts3=NA
  JDW_replustime$xjdsts3=ifelse(JDW_replustime$xjdsts>30,1,ifelse(JDW_replustime$xjdsts<=30,0,"wrong"))
  JDW_replustime$xjdsts3=as.numeric(JDW_replustime$xjdsts3)
  
  save(JDW_replustime,file="JDW_replustime.Rdata")
  write.csv(JDW_replustime,file="JDW_replustime.csv")
  openxlsx::write.xlsx(JDW_replustime,file="JDW_replustime.xlsx")
  # 写并发症
  JDW_replustime_compcalbyima=readxl::read_excel(path ="JDW_replustime11.xlsx",sheet = 1)
  JDW_replustime_compcalbyima_box=JDW_replustime_compcalbyima
  tableshts<- CreateTableOne(vars= c( "coyesno3", "coyesno2"),factorVars =  c("coyesno3","coyesno2"),data = JDW_replustime_compcalbyima_box,includeNA = F, addOverall = T,strata = "xjdsts3")
  tableshts
  dput(names(JDW_replustime_compcalbyima_box))
  c("bah", "name", "diabetes", "complication", "neoadjuvant", "After2020","xjyqts", "xjdsts", "shts", "shts3", "xjyqts3", "xjdsts3", "compcalbyima2","coyesno2", "compcalbyima3", "coyesno3")
  # "coyesno2", 
  JDW_replustime_compcalbyima_box=pivot_wider(
    JDW_replustime_compcalbyima_box,
    id_cols =   c("bah", "name", "diabetes", "complication", "neoadjuvant", "After2020","xjyqts", "xjdsts", "shts", "shts3", "xjyqts3", "xjdsts3"),#coyesno2
    names_from = "compcalbyima2",#compcalbyima3
    names_prefix = "",
    # names_sep = "_",
    names_glue = NULL,
    names_sort = FALSE,
    names_repair = "check_unique",
    values_from = "coyesno2",#coyesno3
    values_fill = NULL,
    values_fn = NULL)
  JDW_replustime_compcalbyima_box$无=NULL
  dput(names(JDW_replustime_compcalbyima_box))
  JDW_replustime_compcalbyima_box[is.na(JDW_replustime_compcalbyima_box)] <- 0
  library(tableone)
  if(F){
  source(r"(C:\Users\zxh\Desktop\R\胃文章\第一篇 全胃wait time与术后关系\发表图表\table\近端胃小篇论文作图表220806\汇总表前的预处理.R)")
  }
  # JDW_replustime_compcalbyima_box=subset(JDW_replustime_compcalbyima_box,bah%in%JDW[which(JDW$pstage1%in%c("III","IV")),"bah"])
  # JDW_replustime_compcalbyima_box=subset(JDW_replustime_compcalbyima_box,bah%in%JDW[which(JDW$After2020==1),"bah"])
  # JDW_replustime_compcalbyima_box=subset(JDW_replustime_compcalbyima_box,bah%in%JDW[which(JDW$neoadjuvant==1),"bah"])
  JDW_replustime_compcalbyima_box=subset(JDW_replustime_compcalbyima_box,bah%in%JDW[which(JDW$neoadjuvant==0),"bah"])
  
  tableshts<- CreateTableOne(vars= c( "吻合口瘘",  "严重反流", "肠梗阻", "切口感染", "腹腔感染", "十二指肠残端瘘", "术后出血", "胸腔积液","胃肠吻合口狭窄"),
    factorVars =  c("吻合口瘘", "严重反流", "肠梗阻", "切口感染", "腹腔感染", "十二指肠残端瘘", "术后出血", "胸腔积液","胃肠吻合口狭窄"),data = JDW_replustime_compcalbyima_box,includeNA = F, addOverall = T,strata = "xjdsts3")
  print(tableshts)
  tablexjdsts=print(tableshts,  formatOptions = list(big.mark = ","),
                   showAllLevels = F, quote = F, noSpaces = T,
                   missing=T, printToggle = T,explain=F)
  write.csv(tablexjdsts,quote=T,file="complication output 2022 in stage III cd2.csv")
  write.csv(tablexjdsts,quote=T,file="complication output 2022 in after 2020 cd2.csv")
  write.csv(tablexjdsts,quote=T,file="complication output 2022 NEO cd2.csv")
  write.csv(tablexjdsts,quote=T,file="complication output 2022 SA cd2.csv")
  write.csv(tablexjdsts,quote=T,file="complication output 2022 cd2.csv")
  
  write.csv(tablexjdsts,quote=T,file="complication output 2022 in stage III.csv",fileEncoding = "UTF-8")
  write.csv(tablexjdsts,quote=T,file="complication output 2022 in after 2020.csv")
  write.csv(tablexjdsts,quote=T,file="complication output 2022 NEO.csv")
  write.csv(tablexjdsts,quote=T,file="complication output 2022 SA.csv")
  write.csv(tablexjdsts,quote=T,file="complication output 2022.csv")

#######################相关性热图  
if(T){ 
    ppp=subset(JDW_replustime,select = c("Sex","Age", "zyts",   "TS", "TSN", "size", "Lauren", "Borrmann", "Diff", "lyvF", "NF", "signet", "pT", "pN", "pM", "eber", "Tumormargin", "surgerymargin","MOR", "NAC","neoadjuvant", "NACT", "NAR", "Smoking", "Drinking", "HD", "hypertension", "diabetes", "Hyperlipidemia", "Antihypertensive", "antidiabetic", "Statins", "sqts", "shts", "PAC", "PACT", "PAR", "rm", "rmt", "death", "dt", "height", "weight", "BMI", "qzryts", "qzssts", "manualsqsjcalculate", "xjyqts","yqts", "xjdsts","dsts","dssj","After2020"))
    ppp_t=subset(JDW_replustime,select = c("neoadjuvant", "NAC","NAR","After2020","xjdsts","xjyqts","zyts","sqts", "shts"))

    pdf(file =paste0("Figure A. 天数之间相关性.pdf"),height = 20,width = 20)
    psych::pairs.panels(ppp_t, 
                        # hist.col="#00FA9A", 
                        show.points=TRUE, 
                        stars=TRUE, 
                        gap=0.05, 
                        pch=".", 
                        ellipses=FALSE, 
                        scale=FALSE,
                        jiggle=TRUE,
                        factor=2,
                        main="Correlation", 
                        # col="#ADFF2F", 
                        pty="m", 
                        cex.labels=3,
                        font=1,
                        # cex.legend=10,
                        cex.cor=1)
    dev.off()
}

  ###
  JDW_replustime_huayan=resultALL_ref_sva[which(resultALL_ref_sva$病案号 %in% JDW_replustime$bah),]
  JDW_replustime_huayan=merge(JDW_replustime_huayan,JDW_replustime,by.x="病案号",by.y = "bah")
  # 限制一下内部的化验
  #如果不写23:59:59 直接写个日期的话右侧不被包含在内
  
  JDW_replustime_huayan=subset(JDW_replustime_huayan,as_date(JDW_replustime_huayan$标本采集时间)>=as_date(JDW_replustime_huayan$szsj)&as_date(JDW_replustime_huayan$标本采集时间)<=as_date(JDW_replustime_huayan$cysj))
  table(as.Date(JDW_replustime_huayan$标本采集时间,tz=Sys.timezone(location = TRUE))<=as.Date(JDW_replustime_huayan$cysj,tz=Sys.timezone(location = TRUE)))
  table(as.Date(JDW_replustime_huayan$标本采集时间,tz=Sys.timezone(location = TRUE))>=as.Date(JDW_replustime_huayan$szsj,tz=Sys.timezone(location = TRUE)))
  length(unique(JDW_replustime_huayan$病案号))
  ###############################################aftersurgery
  JDW_replustime_huayan_aftersurgery=subset(JDW_replustime_huayan,as.Date(JDW_replustime_huayan$标本采集时间,tz=Sys.timezone(location = TRUE))>as.Date(JDW_replustime_huayan$ssrq,tz=Sys.timezone(location = TRUE)))
  length(unique(JDW_replustime_huayan_aftersurgery$病案号))#553
  # 大小
  JDW_replustime_huayan_aftersurgery_nearest=arrange(JDW_replustime_huayan_aftersurgery,desc(JDW_replustime_huayan_aftersurgery$项目结果))
  # JDW_replustime_huayan_aftersurgery_nearest=arrange(JDW_replustime_huayan_aftersurgery,JDW_replustime_huayan_aftersurgery$项目结果)
  JDW_replustime_huayan_aftersurgery_nearest=JDW_replustime_huayan_aftersurgery_nearest[!duplicated(JDW_replustime_huayan_aftersurgery_nearest[,c("项目名称","病案号")]),]
  #POD日期计算
  JDW_replustime_huayan_aftersurgery_nearest$POD=as_date(JDW_replustime_huayan_aftersurgery_nearest$标本采集时间)-JDW_replustime_huayan_aftersurgery_nearest$ssrq
###############
    box=JDW_replustime_huayan_aftersurgery_nearest
    # box$项目结果=log2(box$项目结果+1)
    box=pivot_wider(
      box,
      id_cols = c("病案号","xjdsts3"),#,"zyts","sqts", "shts" ,"qzryts", "qzssts","xjyqts", "xjdsts","yqts", "dsts","neoadjuvant","After2020"
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
    
    select_marker=read.csv(r"(C:\Users\zxh\Desktop\R\jianyanke\用来写术后指标与中文注释.csv)",na.strings = "")
    # box=subset(box,box$项目名称%in% select_marker$项目名称)
if(T){
    wwww1 <- list()
    wwww2 <- list()  
    for(i in 1:5){
    wwww1[[i]]=lapply(subset(res[[i]],select=select_marker$项目名称),function(x) rbindlist(tapply(x,res[[i]]$xjdsts3,shapiro.test))$p.value)
    wwww2[[i]]=lapply(subset(res[[i]],select=select_marker$项目名称),function(x) bartlett.test(x,res[[i]]$xjdsts3)$p.value)
    }
    wwww1=rbindlist(lapply(wwww1,function(x)as.data.frame(t(as.data.frame(x)))))
    wwww2=unlist(wwww2)
    wwww=cbind(wwww1,wwww2)
    wwww$zt=ifelse(wwww$V1>0.05&wwww$V2>0.05,1,0)
    wwww$zt=as.numeric(wwww$zt)
    wwww$fcq=ifelse(wwww$wwww2>0.05,1,0)
    wwww$fcq=as.numeric(wwww$fcq)
    wwww1=NULL
    wwww2=NULL
    # shapiro.test#3-5000
    # lillie.test#>5000
    # tapply(JDW_replustime$zyts,shapiro.test)#不正态
    # bartlett.test(ENSG00000171889~cnv,data=box)#方差齐性
    # shapiro.test(as.numeric(JDW_replustime$dsts))  
    }
    for(i in 1:3){
      library(hrbrthemes)#ggplot2的主题和相关组件包
      library(viridis) #是Matplotlib的新默认颜色映射
      library(scales)

      box=subset(JDW_replustime_huayan_aftersurgery_nearest,JDW_replustime_huayan_aftersurgery_nearest$项目名称 %in%select_marker$项目名称)%>%mutate(`项目名称`=factor(项目名称,levels = select_marker$项目名称))
      box=subset(box,box$病案号%in%rownames(res[[i]]))
      box$xjdsts3=factor(box$xjdsts3,levels = c(0,1))
      box$项目结果=log10(box$项目结果+1)
  p <- ggplot(box, aes(x=xjdsts3, y=`项目结果`)) +
        geom_violin(trim=F,aes(fill=factor(xjdsts3))) +
        #"trim"为TRUE(默认值),将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
        geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.size = 0.001,outlier.color = "yellow",outlier.shape = 24)+ #绘制箱线图
    theme_void()+#+scale_color_manual(values=c("1", "2"))
        # scale_fill_viridis(discrete=TRUE) + #生成一个离散的调色板
        #背景变为白色
        theme(panel.grid.major = element_blank(),   
              panel.grid.minor = element_blank(),
              legend.position="none", #不加图例
              #不显示网格线
              axis.title=  element_blank(),
              panel.spacing = unit(-2, "lines"),
              plot.margin = unit(c(0,0,0,0), "lines"),
              panel.border = element_blank(),
              plot.background = element_blank(),
              # axis.text.y = element_blank(),
              # axis.text.y = element_text(margin=margin(0,0,0,0,"lines")),
              axis.text.x.top  =element_blank(),
              axis.text.x.bottom  = element_text(size=10,color = "black"),#face="bold", #margin=margin(0,0,0,0,"lines"),
              strip.placement = "outside",
              strip.text =element_text(size = 6,angle = 90),#element_blank(),# 
              strip.background = element_blank(),
               # = element_text(,size=2)
              # element_rect(fill = c("red","green")),#
              #去除外框线

              axis.line = element_line(colour = "black",size=0.2)
              #将x轴和y轴加粗显示
        )+labs(title="",x="", y = "log10(Laboratory results + 1)")+stat_compare_means(aes(label = ..p.signif..),method = "kruskal.test", label.y =4,hide.ns =T,size=6,angle=90)+facet_grid(rows = vars(box$项目名称),margins = F,switch = "both")+coord_flip(xlim = c(0,4.2)) #翻转坐标
      
    pdf(file =paste0("Figure A.",i,"violin.pdf"),height = 10,width = 1.3)
    print(p)
    dev.off()
    openxlsx::write.xlsx(ggplot_build(p)$data,file =paste0("Figure A.",i,"violin.xlsx"))
}
    ####################
    
    
    
    
    
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
    pp11=pp11[which(pp11$xName%in%c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "xjyqts", "xjdsts")|pp11$yName%in%c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "xjyqts", "xjdsts")),]
    pp11=subset(pp11,xName!=yName)
    pp11=subset(pp11,!pp11$yName%in%c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "xjyqts", "xjdsts"))
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
    ts_husyan_sign22=subset(ts_husyan_sign,!ts_husyan_sign$yName%in%c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "xjyqts", "xjdsts"))
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
  ts_husyan_sign33=subset(ts_husyan_sign,!ts_husyan_sign$yName%in%c("zyts","sqts", "shts" ,"qzryts", "qzssts",  "xjyqts", "xjdsts"))
  ts_husyan_sign33=subset(ts_husyan_sign33,abs(sign(ts_husyan_sign33$corr)+sign(ts_husyan_sign33$corr.1))==2)
  # ts_husyan_sign33=subset(ts_husyan_sign33,abs(ts_husyan_sign33$corr)>0.2&abs(ts_husyan_sign33$corr.1)>0.2)
  openxlsx::write.xlsx(ts_husyan_sign33,file ="ts_husyan_sign33.xlsx")
  ######################
  
  # ts_husyan_sign22=subset(ts_husyan_sign22,abs(sign(ts_husyan_sign22$corr)+sign(ts_husyan_sign22$corr.1)+sign(ts_husyan_sign22$corr.2))==3)

  box=JDW_replustime_huayan_aftersurgery_nearest
  # box$项目结果=log2(box$项目结果+1)
  box=tidyr::pivot_longer(
      data=box,
      cols=c("xjyqts", "xjdsts","zyts", "shts"),#everthing,"sqts","qzryts", "qzssts",  
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
    pdf(file =paste0("Figure A. 化验全因相关.pdf"),height = 200,width = 20)
    box$Daytime=as.numeric(box$Daytime)
    # box$Daytime=log2(box$Daytime+1)
    p1=ggplot(box, aes(x = 项目结果, y = Daytime)) + 
      ylab("")+xlab("")+
      geom_point(shape = 21, colour = "#4682B4", fill = "#87CEFA", size = 3, stroke = .5,alpha=0.8)+ geom_smooth(method="glm",formula = y ~ x,linetype=2,color="#6495ED",fill="#D3D3D3") + theme_bw()+stat_cor(method = 'spearman', aes(),na.rm = T,size = 10)+facet_grid(rows=vars(box$项目名称),cols=vars(box$DAYtype),margins = F)#,space="free_y"
    # p2=ggExtra::ggMarginal(p1, type = "density", xparams = list(fill = "#FFE4B5"),yparams = list(fill = "#90EE90"))
    p1
    dev.off()



  box=JDW_replustime_huayan_aftersurgery_nearest
  write.csv(JDW_replustime_huayan_aftersurgery_nearest,file="box.csv")
  # box$项目结果=log2(box$项目结果+1)
  box=pivot_wider(
    box,
    id_cols = c("病案号","neoadjuvant","Sex", "Age", "TS", "size", "Lauren", "Borrmann", "Diff", "signet","pT","pN","pM","Smoking", "Drinking", "HD", "hypertension", "diabetes","Age","Sex","pT","ssrq","zyts","sqts", "shts" ,"qzryts", "qzssts","yqts", "xjdsts","zyts1","sqts1", "shts1" ,"qzryts1", "qzssts1",  "xjyqts1", "xjyqts2", "xjyqts3", "xjdsts1","xjdsts2","xjdsts3"),#everthing#,"neoadjuvant","After2020"
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
  
  #开始建模
  JDW_rrrrr=subset(JDW_replustime,select=c("bah","Sex", "Age","TS", "size", "Lauren", "Borrmann", "Diff", "lyvF", "NF", "signet", "pT", "pN",  "pM", "Tumormargin", "surgerymargin", "Surgicalapproach", "Surgicalway", "NAC", "NACT", "NAR", "Smoking", "Drinking", "HD", "hypertension", "diabetes", "Hyperlipidemia",   "Antihypertensive", "antidiabetic", "Statins", "shts3", "PAC", "PAR","height", "weight", "BMI","neoadjuvant", "After2020","xjdsts3","xjdsts"))
  JDW_rrrrr$BMI=JDW_rrrrr$weight/(JDW_rrrrr$height/100)^2
  JDW=JDW_rrrrr
  JDW$Age=as.integer(JDW$Age)
  JDW$Age2 = cut(JDW$Age,breaks = c(-Inf,65,Inf),right =T,ordered_result = T,
                 labels = c("<65","≥65"))
  JDW$TS=ifelse(is.na(JDW$TS),NA,
                ifelse(JDW$TS%in%c("0","1","2"),"proximal",
                       ifelse(JDW$TS%in%c("3","4","5"),"middle/distal",NA)))
  JDW$Lauren= ifelse(is.na(JDW$Lauren),NA,
                     ifelse(JDW$Lauren=="0","0",
                            ifelse(JDW$Lauren=="1","1",
                                   ifelse(JDW$Lauren=="2","2",NA))))
  JDW$Lauren=factor(JDW$Lauren,levels = c(0:2),ordered = F,
                    labels=c("intestinal","mixed","diffuse"))
  JDW[JDW$Borrmann=="5","Borrmann"]=NA
  JDW$Borrmann=as.integer(  JDW$Borrmann)
  # JDW$Borrmann=ifelse(is.na(JDW$Borrmann),NA,
  #                     ifelse(JDW$Borrmann%in%c("0","1"),"superficial",
  #                            ifelse(JDW$Borrmann%in%c("2","3","4"),"ulcerative",NA)))
  JDW$Diff=factor(JDW$Diff,
                  levels = c(0:5),
                  labels = c("Poorly differentiated","Poorly differentiated",
                             "Poorly differentiated","Well differentiated",
                             "Well differentiated","Well differentiated"),ordered = T)
  
  
  # 接着转变数据类型和清洗数据
  # JDW$lyvF=as.character(JDW$lyvF)
  # JDW$NF=as.character(JDW$NF)
  # JDW$signet=as.character(JDW$signet)
  # JDW$lyvF=ifelse(is.na(JDW$lyvF),NA,
  #                 ifelse(JDW$lyvF=="0","Negative",
  #                        ifelse(JDW$lyvF=="1","Positive",NA)))
  # JDW$NF=ifelse(is.na(JDW$NF),NA,
  #               ifelse(JDW$NF=="0","Negative",
  #                      ifelse(JDW$NF=="1","Positive",NA)))
  # JDW$signet=factor(JDW$signet,levels = c(0:2),ordered = F,
  #                   labels = c("No-Signet","Partial-Signet","All-Signet"))
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
  JDW$pstage = ifelse(
    JDW$pM == "1","IV",
    ifelse(JDW$tplusn =="1", "IA",
           ifelse(JDW$tplusn =="2", "IB",
                  ifelse(JDW$tplusn =="3", "IIA",
                         ifelse(JDW$tplusn =="4", "IIB",
                                ifelse(JDW$tplusn =="5", "IIIA",
                                       ifelse(JDW$tplusn=="6" & JDW$pT1=="4" , "IIIA",
                                              ifelse(JDW$tplusn%in%c("6","7"), "IIIB", 
                                                     ifelse(JDW$tplusn%in%c("8","9","10"), "IIIC", NA)))))))))
  table(JDW$tplusn)
  table(JDW$pstage)
  JDW$pstage1=ifelse(is.na(JDW$pstage),NA,
                     ifelse(JDW$pstage %in%c("IA", "IB"),"I",
                            ifelse(JDW$pstage %in%c("IIA","IIB"), "II",
                                   ifelse(JDW$pstage %in%c("IIIA","IIIB","IIIC"), "III",
                                          ifelse(JDW$pstage =="IV", "IV",NA)))))
  
  
  JDW$pT=ifelse(is.na(JDW$pT),NA,
                ifelse(JDW$pT%in%c("0","1","1a","1b"),"early",
                       ifelse(JDW$pT%in%c("2","2a","2b",
                                          "3","3a","3b","4","4a","4b"),"advanced",NA)))
  JDW$pN=ifelse(is.na(JDW$pN),NA,
                ifelse(JDW$pN=="0","Negative",
                       ifelse(JDW$pN%in%c("1","1a","1b","2",
                                          "2a","2b","3","3a","3b"),"Positive",NA)))
  JDW$pM01=as.numeric(JDW$pM)
  JDW$pM=ifelse(is.na(JDW$pM),NA,
                ifelse(JDW$pM=="0","Negative",
                       ifelse(JDW$pM%in%c("1"),"Positive",NA)))
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
  # JDW$sqts=as.integer(JDW$sqts)
  # JDW$shts=as.integer(JDW$shts)
  # JDW$yqts=as.integer(JDW$yqts)
  # JDW$qzryts=as.integer(JDW$qzryts)
  # JDW$qzssts=as.integer(JDW$qzssts)
  # JDW$dsts=as.integer(JDW$dsts)
  # JDW$NAC=as.character(JDW$NAC)
  # JDW$NAR=as.character(JDW$NAR)
  # JDW$PAC=as.character(JDW$PAC)
  # JDW$PAR=as.character(JDW$PAR)
  dput(names(JDW)) # 输出据集变量名称
  str(JDW)
  JDW$pT=factor(JDW$pT,levels = c("early","advanced"),ordered = T)
  JDW$Sex=ifelse(JDW$Sex==1,"Male","Female")
  JDW$xjdsts3=factor(JDW$xjdsts3,levels = c(0,1),labels=c("Normal","Delayed"),ordered = T)
  JDW <- Units(JDW,list("Age"="odds per year"))
  uni.odds <- glmSeries(shts3~1,vars=c( "Sex", "Age","TS", "size", "Lauren", "Borrmann", "Diff", "lyvF", "NF", "signet", "pT", "pN",  "pM01", "surgerymargin", "Surgicalapproach", "Surgicalway", "NAC", "NACT", "NAR", "Smoking", "Drinking", "HD", "hypertension", "diabetes", "Hyperlipidemia", "pstage","pstage1",  "Antihypertensive", "antidiabetic", "Statins", "PAC", "PAR","height", "weight", "BMI","neoadjuvant", "After2020","xjdsts3"),data=JDW,family=binomial())
  uni.odds <- glmSeries(shts3~1,vars=c("xjdsts3","Sex", "Age","TS", "pT", "pN","Smoking"),data=JDW,family=binomial())
 publish(uni.odds)
  write.csv(uni.odds,file = "uni.odds.csv")
  nrow(na.omit(subset(JDW,select=c("xjdsts3","Sex", "Age","TS", "pT", "pN","Smoking"))))#584
  x=uni.odds
  num <- length(names(x))
  plotcols<-x[,(num-3):(num-1)]
  tabcols <-subset(x,select=c(1,2,6))
  tabcols$Pvalue=round(tabcols$Pvalue,2)
  eee=Publish::plotConfidence(x=plotcols, labels=tabcols)
  uni.odds=uni.odds%>%mutate(across(where(is.numeric),~round(.,2)))

  write.csv(  cbind(uni.odds,eee[["values"]][["labels"]]),file = "uni.odds.csv")
  
  
  f1 <- glm(shts3~xjdsts3+Sex+Age+TS+pT+pN+Smoking,data=JDW,family=binomial())
  publish(f1)
  plot(regressionTable(f1,factor.reference = "inline"),cex=1.3)
  summary(regressionTable(f1,factor.reference = "extraline"))$regressionTable
  eee=summary(regressionTable(f1,factor.reference = "inline"))$regressionTable
  eee=tidyr::unite(eee,col="Multivariable Adjusted Model",OddsRatio,CI.95,remove = TRUE,sep = "")
  write.csv(eee,file = "multiple.csv")


  # summary(regressionTable(fit),handler="prettyNum")
  # summary(regressionTable(fit),handler="format")
  # summary(regressionTable(fit),handler="sprintf",digits=c(2,2),pValue.stars=TRUE)
  # summary(regressionTable(fit),handler="sprintf",digits=c(2,2),pValue.stars=TRUE,ci.format="(l,u)")
  
  # JDW_replustime_compcalbyima$coyesno3#有无
  # JDW_replustime_compcalbyima$compcalbyima3#具体
  # 没有意义
  if(F){

    JDW_replustime_compcalbyima_box_log=merge(subset(JDW_replustime_compcalbyima_box,select=c(bah,coyesno3)),JDW,by="bah")
  coyesno3.odds <- glmSeries(coyesno3~1,vars=c( "Sex", "Age","TS", "size", "Lauren", "Borrmann", "Diff", "lyvF", "NF", "signet", "pT", "pN",  "pM01", "surgerymargin", "Surgicalapproach", "Surgicalway", "NAC", "NACT", "NAR", "Smoking", "Drinking", "HD", "hypertension", "diabetes", "Hyperlipidemia", "pstage","pstage1",  "Antihypertensive", "antidiabetic", "Statins", "PAC", "PAR","height", "weight", "BMI","neoadjuvant", "After2020","xjdsts3"),data=JDW_replustime_compcalbyima_box_log,family=binomial())
  coyesno3.odds
  coyesno3.odds <- glmSeries(JDW_replustime_compcalbyima_box_log~1,vars=c("xjdsts3","Sex", "Age","TS", "pT", "pN","Smoking"),data=JDW,family=binomial())
}
  ## control the logistic regression analyses for age and gender 
  ## but collect only information on the variables in `vars'.

  if(F){
    fit_p <- glm(xjdsts3~Age+Sex+Lauren,family=binomial(),data=JDW_replustime)
    summary(fit_p)#效果不佳
    publish(fit_p)#好
    rtab <- Publish::regressionTable(fit_p,factor.reference = "extraline")#inline
    summary(rtab)#好
    
    
    controlled.odds <- glmSeries(hyper1~Age+Sex+TS+,
                                 vars=c("chol","hdl","location"),
                                 data=JDW_replustime, family=binomial)
    controlled.odds
    JDW_replustime
    
    
  }
  
  # publish(t.test(bp.2s~gender,data=Diabetes))
  # publish(wilcox.test(bp.2s~gender,data=Diabetes))
  # publish(with(Diabetes,t.test(bp.2s,bp.1s,paired=TRUE)))
  # publish(with(Diabetes,wilcox.test(bp.2s,bp.1s,paired=TRUE)))
  plot(rtab,cex=1.8)
  plot(rtab,xlim=c(0.85,1.15),cex=1.8,xaxis.cex=1.5)
  fit_p <- glm(xjdsts3~GGT+AST+TT+MCV,family=binomial(),data=box)
  summary(fit_p)
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
  ggplot(JDW, aes(x = xjyqts, fill = NAC)) +
    # 直方图函数：position设置堆积模式为重叠
    geom_histogram(position = "identity", alpha = 0.4)
  table(JDW_replustime_huayan_nearest$xjyqts,JDW_replustime_huayan_nearest$NAC)
  boxplot(xjyqts~NAC,JDW)
  
  

  
  
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
  