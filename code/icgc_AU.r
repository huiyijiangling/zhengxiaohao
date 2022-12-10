#DFI: 
#1 
#for patient having new tumor event whether it is a local recurrence, distant metastasis, new primary tumor of the cancer, including cases with a new tumor event whose type is N/A.  
#Disease free was defined by: first, treatment_outcome_first_course is "Complete Remission/Response";
#if the tumor type doesn't have "treatment_outcome_first_course" then disease-free was defined by the value "R0" in the field of "residual_tumor"; 
#otherwise, disease-free was defined by the value "negative" in the field of "margin_status". 
#If the tumor type did not have any of these fields, then its DFI was NA. 指的是切除后完整性
#0
#New primary tumor in other organ was censored;
#patients who were Dead with tumor without new tumor event are excluded;
#patients wih stage IV are excluded too.
#analysis
# os >=3month,


rm(list=ls())
options(stringsAsFactors = F)
gc()
# ICGC_PACA_AU = read_tsv(file="C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/ICGC/guanwang/PACA-AU/exp_seq.tsv")
# View(ICGC_PACA_AU)ICGC_PACA_AU
# system("cat C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/ICGC/guanwang/PACA-AU/exp_seq.tsv.gz |cut -f 1,3,8,9 >exp_seq_PRAD-CA_simplify.txt")
library(readr)
library(dplyr)
library(tidyverse)

donor= read_tsv(file="./ICGC/guanwang/PACA-AU/donor.tsv.gz",guess_max=min(1000000, Inf))#read_delim
specimen= read_tsv(file="./ICGC/guanwang/PACA-AU/specimen.tsv.gz",guess_max=min(1000000, Inf))#read_delim
sample=read_tsv(file="./ICGC/guanwang/PACA-AU/sample.tsv.gz",guess_max=min(1000000, Inf))#read_delim
parsing_failures <- problems(specimen)
parsing_failures
# kkkkkk<- list(donor, sample, specimen) %>% reduce(full_join, by = "icgc_donor_id") #这样其实产生了很多错误
# kkkkkk<- list(sample, specimen) %>% reduce(full_join, by = "icgc_specimen_id")
sample_specimen<-merge(sample,specimen,by="icgc_specimen_id",all=T)
sample_specimen_donor<-merge(sample_specimen,donor,by.x="icgc_donor_id.x",by.y="icgc_donor_id",all=T)
ICGC_PACA_AU_metaMatrix=sample_specimen_donor
table(duplicated(sample_specimen_donor$icgc_specimen_id))
table(duplicated(sample_specimen_donor$icgc_donor_id.x))
table(duplicated(sample_specimen_donor$icgc_sample_id))
# kkkkkk2<- list(sample, donor, specimen) %>% reduce(left_join, by = "icgc_donor_id")

#处理seq
ICGC_PACA_AU_seq = read_tsv(file="./ICGC/guanwang/PACA-AU/exp_seq.tsv.gz",guess_max=min(1000000, Inf))#read_delim
## 查看sample 数目
length(unique(ICGC_PACA_AU_seq$icgc_donor_id))    # [1] 91
length(unique(ICGC_PACA_AU_seq$icgc_specimen_id)) # [1] 92
length(unique(ICGC_PACA_AU_seq$icgc_sample_id)) # [1] 92 就用他
#icgc_sample_id是唯一的
#icgc_specimen_id可能重复
#icgc_donor_id是重复的
#spread只能留一个列，多了否则就会报错
# gather 和spread 就是结构和合并的方法
ICGC_PACA_AU_seq_norm=ICGC_PACA_AU_seq[,c("normalized_read_count","gene_id","icgc_sample_id")]
ICGC_PACA_AU_seq_raw=ICGC_PACA_AU_seq[,c("raw_read_count","gene_id","icgc_sample_id")]
ICGC_PACA_AU_seq_norm_test=unite(data=ICGC_PACA_AU_seq_norm,col = "e",c("gene_id","icgc_sample_id"),remove = F)
table(duplicated(ICGC_PACA_AU_seq_norm_test$e))#
#
ICGC_PACA_AU_seq_norm <- ICGC_PACA_AU_seq_norm %>% group_by(gene_id)  %>% spread(icgc_sample_id,normalized_read_count) 
ICGC_PACA_AU_seq_norm=as.data.frame(ICGC_PACA_AU_seq_norm)
rownames(ICGC_PACA_AU_seq_norm)=ICGC_PACA_AU_seq_norm$gene_id
ICGC_PACA_AU_seq_norm=select(ICGC_PACA_AU_seq_norm,-"gene_id")
ICGC_PACA_AU_seq_raw <- ICGC_PACA_AU_seq_raw %>% group_by(gene_id) %>% spread(icgc_sample_id,raw_read_count)
ICGC_PACA_AU_seq_raw=as.data.frame(ICGC_PACA_AU_seq_raw)
rownames(ICGC_PACA_AU_seq_raw)=ICGC_PACA_AU_seq_raw$gene_id
ICGC_PACA_AU_seq_raw=select(ICGC_PACA_AU_seq_raw,-"gene_id")
#####################################
library(edgeR)
expr = DGEList(counts = ICGC_PACA_AU_seq_raw)
keepALL <- rowSums(cpm(expr) > 1) >= 0.5*ncol(ICGC_PACA_AU_seq_raw)
nGenes <- as.numeric(summary(keepALL)[2]) + 
  as.numeric(summary(keepALL)[3])
nKeep <- summary(keepALL)[3]
cat (paste('Total Number of genes: ', nGenes, '\n', sep=''))
cat (paste('Number of genes for downstream analysis: ', nKeep, 
           '\n', sep=''))
exprSet=ICGC_PACA_AU_seq_raw[keepALL,]
# exprSet=rnaCounts_Pancreas[,group_list=="tumor"][keepALL,]
# dds <- dds[keep, ]
# dge <- dge[keep,,keep.lib.sizes = TRUE]
exprSet_rownames <- rownames(exprSet)
exprSet_colnames <- colnames(exprSet)
exprSet_quant <- preprocessCore::normalize.quantiles(
  as.matrix(exprSet))
rownames(exprSet_quant) <- exprSet_rownames
colnames(exprSet_quant) <- exprSet_colnames
ICGC_PACA_AU_seq_log_filter=log2(exprSet_quant+1)
if(T){
#处理array
ICGC_PACA_AU_array = read_tsv(file="./ICGC/guanwang/PACA-AU/exp_array.PACA-AU.tsv.gz")#read_delim
## 查看sample 数目
length(unique(ICGC_PACA_AU_array$icgc_donor_id))    # [1] 91
length(unique(ICGC_PACA_AU_array$icgc_specimen_id)) # [1] 92
length(unique(ICGC_PACA_AU_array$icgc_sample_id)) # [1] 92 就用他
#icgc_sample_id,icgc_specimen_id是唯一的
#icgc_donor_id是重复的
#spread只能留一个列，多了否则就会报错
# gather 和spread 就是结构和合并的方法
ICGC_PACA_AU_array_norm=ICGC_PACA_AU_array[,c("normalized_expression_value","gene_id","icgc_sample_id")]
#
ICGC_PACA_AU_array_norm <- ICGC_PACA_AU_array_norm %>% group_by(gene_id)  %>% spread(icgc_sample_id,normalized_expression_value)
ICGC_PACA_AU_array_norm=as.data.frame(ICGC_PACA_AU_array_norm)
rownames(ICGC_PACA_AU_array_norm)=ICGC_PACA_AU_array_norm$gene_id
ICGC_PACA_AU_array_norm=select(ICGC_PACA_AU_array_norm,-"gene_id")#GSE36924 可以还原

library(AnnoProbe)
(gpl="GPL10558")
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl=gpl,type = "bioc")#soft的合并性差落后，建议首选bioc,但如果结合了更新就不一样了,例如ta
head(probe2gene)
source("./updateName.R")# ncbi与biomart相比多map了200，等于是，ncbi》biomart》org.hs
forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol)
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
probe2gene=unique(probe2gene[,c("probe_id","ENSEMBL")])#bioc
# probe2gene=unique(probe2gene[,c("ID","ENSEMBL")])#soft
# probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
colnames(probe2gene)=c("probe_id","symbol")
ICGC_PACA_AU_array_norm <- filterEM(ICGC_PACA_AU_array_norm,probe2gene)
}
write.csv(ICGC_PACA_AU_metaMatrix,file="ICGC_PACA_AU_metaMatrix.csv",quote = T)
save(ICGC_PACA_AU_seq_log_filter,ICGC_PACA_AU_seq_norm,ICGC_PACA_AU_seq_raw,ICGC_PACA_AU_seq,ICGC_PACA_AU_array_norm,ICGC_PACA_AU_metaMatrix,file = "ICGC_PACA_AU_expr.Rdata")
load("ICGC_PACA_AU_expr.Rdata")
#去掉杂质
#去掉经过术前治疗病历 保留 无治疗 手术 na 
table(ICGC_PACA_AU_metaMatrix$specimen_donor_treatment_type,useNA = "ifany")
ICGC_PACA_AU_metaMatrix=subset(ICGC_PACA_AU_metaMatrix,specimen_donor_treatment_type %in% c("surgery","no treatment"))#不能包含NA，否则有可能进来化疗后的
#duct
table(ICGC_PACA_AU_metaMatrix$specimen_type,useNA = "ifany")
table(ICGC_PACA_AU_metaMatrix$specimen_type_other,useNA = "ifany")
table(ICGC_PACA_AU_metaMatrix$tumour_histological_type,useNA = "ifany")
ICGC_PACA_AU_metaMatrix_ductNA=subset(ICGC_PACA_AU_metaMatrix,tumour_histological_type %in% c(NA,"Pancreatic Ductal Adenocarcinoma"))
ICGC_PACA_AU_metaMatrix_ductNA_TAN=subset(ICGC_PACA_AU_metaMatrix_ductNA,specimen_type %in% 
                                            c("Normal - solid tissue","Normal - tissue adjacent to primary","Primary tumour - solid tissue","Primary tumour - other"))#不能房间来NA，不是就不是
# Primary tumour - other 有些肿瘤可能很重要
ICGC_PACA_AU_metaMatrix_ductNA_T=subset(ICGC_PACA_AU_metaMatrix_ductNA_TAN,specimen_type %in% c("Primary tumour - solid tissue"))
ICGC_PACA_AU_metaMatrix_ductNA_A=subset(ICGC_PACA_AU_metaMatrix_ductNA_TAN,specimen_type %in% c("Normal - tissue adjacent to primary"))
ICGC_PACA_AU_metaMatrix_ductNA_N=subset(ICGC_PACA_AU_metaMatrix_ductNA_TAN,specimen_type %in% c("Normal - solid tissue"))
table(ICGC_PACA_AU_metaMatrix_ductNA_T$tumour_histological_type,useNA = "ifany")
table(ICGC_PACA_AU_metaMatrix_ductNA_A$tumour_histological_type,useNA = "ifany")
table(ICGC_PACA_AU_metaMatrix_ductNA_N$tumour_histological_type,useNA = "ifany")
#
ICGC_PACA_AU_metaMatrix_duct_T=subset(ICGC_PACA_AU_metaMatrix_ductNA_T,tumour_histological_type %in% c("Pancreatic Ductal Adenocarcinoma"))
table(ICGC_PACA_AU_metaMatrix_duct_T$tumour_histological_type,useNA = "ifany")
#
# table(ICGC_PACA_AU_metaMatrix_duct_T$donor_vital_status)                   # counts the "unusual" NA 等于  table(ICGC_PACA_AU_metaMatrix_duct_T$donor_vital_status, exclude = NA)     # counts none
# table(ICGC_PACA_AU_metaMatrix_duct_T$donor_vital_status,useNA = "ifany")  # counts all three 等于  table(ICGC_PACA_AU_metaMatrix_duct_T$donor_vital_status, exclude = NULL)   #  (ditto)
table(ICGC_PACA_AU_metaMatrix_duct_T$tumour_histological_type,useNA = "ifany")
# OS
ICGC_PACA_AU_metaMatrix_duct_T$OS.time=apply(ICGC_PACA_AU_metaMatrix_duct_T[,c("donor_relapse_interval","donor_interval_of_last_followup","donor_survival_time")],1,max)
ICGC_PACA_AU_metaMatrix_duct_T$OSbi=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$donor_vital_status),NA,
                                                ifelse(ICGC_PACA_AU_metaMatrix_duct_T$donor_vital_status=="alive","0",
                                                       ifelse(ICGC_PACA_AU_metaMatrix_duct_T$donor_vital_status=="deceased","1","WRONG")))
ICGC_PACA_AU_metaMatrix_duct_T$OS=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$OSbi)&is.na(ICGC_PACA_AU_metaMatrix_duct_T$OS.time),NA,
                                         ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$OSbi)&(!is.na(ICGC_PACA_AU_metaMatrix_duct_T$OS.time)),"0",
                                                ifelse(!is.na(ICGC_PACA_AU_metaMatrix_duct_T$OSbi),ICGC_PACA_AU_metaMatrix_duct_T$OSbi,"WRONG")))
# 某些人登记错的，rfi比OS还大
table(ICGC_PACA_AU_metaMatrix_duct_T$donor_relapse_interval>ICGC_PACA_AU_metaMatrix_duct_T$OS.time)
# Relapse
# official tools has many problem 直接把relpase，progression当做复发，其他的sd和nes，和pr均当做不复发
#注意xx$xx时用==，na是不被考虑的
table(ICGC_PACA_AU_metaMatrix_duct_T$donor_relapse_interval,useNA = "ifany")
#RFI不需要分RFIbi，我算都一样
ICGC_PACA_AU_metaMatrix_duct_T$RFI=ifelse(!is.na(ICGC_PACA_AU_metaMatrix_duct_T$donor_relapse_interval),"1",
                                          ifelse(!is.na(ICGC_PACA_AU_metaMatrix_duct_T$donor_relapse_type),"1",
                                                 ifelse(ICGC_PACA_AU_metaMatrix_duct_T$disease_status_last_followup %in% c("relapse",'progression','partial remission'),"1",
                                                        ifelse(ICGC_PACA_AU_metaMatrix_duct_T$disease_status_last_followup %in% c("complete remission",'no evidence of disease','stable'),"0",
                                                               ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$disease_status_last_followup),NA,"Wrong")))))
ICGC_PACA_AU_metaMatrix_duct_T$RFI.time=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$RFI),NA,
                                               ifelse(ICGC_PACA_AU_metaMatrix_duct_T$RFI=="0",ICGC_PACA_AU_metaMatrix_duct_T$OS.time,
                                                      ifelse(ICGC_PACA_AU_metaMatrix_duct_T$RFI=="1",ICGC_PACA_AU_metaMatrix_duct_T$donor_relapse_interval,"WRONG")))
ICGC_PACA_AU_metaMatrix_duct_T$DFS=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$RFI),NA,
                                          ifelse(ICGC_PACA_AU_metaMatrix_duct_T$RFI=="0"&ICGC_PACA_AU_metaMatrix_duct_T$OS=="1","1",
                                                 ifelse(!is.na(ICGC_PACA_AU_metaMatrix_duct_T$RFI),ICGC_PACA_AU_metaMatrix_duct_T$RFI,"Wrong")))
ICGC_PACA_AU_metaMatrix_duct_T$DFS.time=ICGC_PACA_AU_metaMatrix_duct_T$RFI.time
table(ICGC_PACA_AU_metaMatrix_duct_T$RFI=="0"&ICGC_PACA_AU_metaMatrix_duct_T$OS=="1")
table(ICGC_PACA_AU_metaMatrix_duct_T$RFI,useNA = "ifany")
table(ICGC_PACA_AU_metaMatrix_duct_T$DFS,useNA = "ifany")
table(ICGC_PACA_AU_metaMatrix_duct_T$DFS,useNA = "ifany")
table(ICGC_PACA_AU_metaMatrix_duct_T$DFS.time,useNA = "ifany")# Relapse type
table(ICGC_PACA_AU_metaMatrix_duct_T$donor_relapse_type,useNA = "ifany")# Relapse type

ICGC_PACA_AU_metaMatrix_duct_T_seestage=distinct(ICGC_PACA_AU_metaMatrix_duct_T[is.na(ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage),c("icgc_sample_id","tumour_stage","tumour_stage_supplemental","donor_tumour_stage_at_diagnosis_supplemental","donor_tumour_stage_at_diagnosis")])
ICGC_PACA_AU_metaMatrix_duct_T$TNM=ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage
# ICGC_PACA_AU_metaMatrix_duct_T[ICGC_PACA_AU_metaMatrix_duct_T$icgc_sample_id %in% "SA533649",]$TNM <- ICGC_PACA_AU_metaMatrix_duct_T[ICGC_PACA_AU_metaMatrix_duct_T$icgc_sample_id %in% "SA533649",]$donor_tumour_stage_at_diagnosis_supplemental
ICGC_PACA_AU_metaMatrix_duct_T$stage123=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$donor_tumour_stage_at_diagnosis_supplemental)& is.na(ICGC_PACA_AU_metaMatrix_duct_T$donor_tumour_stage_at_diagnosis)&is.na(ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage_supplemental)&is.na(ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage),NA,
                                               ifelse(grepl("m1",ICGC_PACA_AU_metaMatrix_duct_T$TNM,ignore.case=T),"no",
                                                      ifelse(grepl("iv",ICGC_PACA_AU_metaMatrix_duct_T$TNM,ignore.case=T),"no","yes")))
# DFI 去掉m1，去掉随访不足1月
table(ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,useNA = "ifany")# Relapse type
table(grepl("m1",ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,ignore.case=T))

if(F){
  #干掉m1
ICGC_PACA_AU_metaMatrix_duct_T$DFS.time=ifelse(ICGC_PACA_AU_metaMatrix_duct_T$stage123=="no",NA,as.numeric(ICGC_PACA_AU_metaMatrix_duct_T$DFS.time))
ICGC_PACA_AU_metaMatrix_duct_T$DFS=ifelse(ICGC_PACA_AU_metaMatrix_duct_T$stage123=="no",NA,ICGC_PACA_AU_metaMatrix_duct_T$DFS)

ICGC_PACA_AU_metaMatrix_duct_T$RFI.time=ifelse(ICGC_PACA_AU_metaMatrix_duct_T$stage123=="no",NA,as.numeric(ICGC_PACA_AU_metaMatrix_duct_T$RFI.time))
ICGC_PACA_AU_metaMatrix_duct_T$RFI=ifelse(ICGC_PACA_AU_metaMatrix_duct_T$stage123=="no",NA,ICGC_PACA_AU_metaMatrix_duct_T$RFI)
}


if(F){
# 理论上还要干掉90天
ICGC_PACA_AU_metaMatrix_duct_T$DFI=ifelse(ICGC_PACA_AU_metaMatrix_duct_T$OS.time>=90,ICGC_PACA_AU_metaMatrix_duct_T$DFS,NA)
ICGC_PACA_AU_metaMatrix_duct_T$DFI.time=ifelse(ICGC_PACA_AU_metaMatrix_duct_T$OS.time>=90,ICGC_PACA_AU_metaMatrix_duct_T$DFS.time,NA)
}
# 理论上还要干掉r0
#不建议放这里处理
# ghjkl=subset(ICGC_PACA_AU_metaMatrix_duct_T,grepl("m1",ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,ignore.case=T))
ICGC_PACA_AU_metaMatrix_duct_surv=ICGC_PACA_AU_metaMatrix_duct_T[,c("icgc_donor_id.x","DFS","DFS.time","OS","OS.time","RFI","RFI.time","stage123")]#,"tumour_stage"
ICGC_PACA_AU_metaMatrix_duct_surv=dplyr::distinct(ICGC_PACA_AU_metaMatrix_duct_surv)
# tapply(ppp$OS.time,ppp$specimen_donor_treatment_type,median)
# qqq=ICGC_PACA_AU_metaMatrix_duct_T[grep("mx",ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,ignore.case=T),]

library(AnnoProbe)
probe2gene=unique(ICGC_PACA_AU_metaMatrix_duct_T[,c("icgc_sample_id","icgc_donor_id.x")])#soft
colnames(probe2gene)=c("probe_id","symbol")
sampletodonor <- function(probes_expr){
  probes_expr=as.data.frame(t(probes_expr))
  genes_expr <- filterEM(probes_expr,probe2gene)
  genes_expr=as.data.frame(t(genes_expr))
  return(genes_expr)
}
ICGC_PACA_AU_seq_norm_dornor=sampletodonor(ICGC_PACA_AU_seq_norm)#92-58
ICGC_PACA_AU_seq_raw_dornor=sampletodonor(ICGC_PACA_AU_seq_raw)#92-58
ICGC_PACA_AU_array_norm_dornor=sampletodonor(ICGC_PACA_AU_array_norm)#269-230
ICGC_PACA_AU_seq_log_filter_dornor=sampletodonor(ICGC_PACA_AU_seq_log_filter)#92-58

save(ICGC_PACA_AU_seq_log_filter_dornor,ICGC_PACA_AU_seq_norm_dornor,ICGC_PACA_AU_seq_raw_dornor,ICGC_PACA_AU_array_norm_dornor,ICGC_PACA_AU_metaMatrix_duct_surv,ICGC_PACA_AU_metaMatrix_duct_T,file = "ICGC_PACA_AU_metaMatrix_duct_surv.Rdata")
# ICGC_PACA_AU_seq_log_filter_dornor 如果是建立risk score 没必要因为实际上array里的人是重复的
load("ICGC_PACA_AU_metaMatrix_duct_surv.Rdata")
write.csv(ICGC_PACA_AU_metaMatrix_duct_surv,file="ICGC_PACA_AU_metaMatrix_duct_surv.csv",quote = T)
write.csv(ICGC_PACA_AU_metaMatrix_duct_T,file="ICGC_PACA_AU_metaMatrix_T.csv",quote = T)

#是不是应该 keepsize


# ICGC_PACA_AU_metaMatrix_duct_surv
# 
# 
# 
# DO32932
# 
# 
# 
# 
# 
# v5[duplicated(v5$icgc_donor_id.x),"icgc_donor_id.x"]
# 
# 
# jjjjjj=ICGC_PACA_AU_metaMatrix_duct_surv[order(ICGC_PACA_AU_metaMatrix_duct_surv$icgc_donor_id.x,decreasing = T),]
# 
# gg=ICGC_PACA_AU_metaMatrix_duct_surv[duplicated(ICGC_PACA_AU_metaMatrix_duct_surv$icgc_donor_id.x),"icgc_donor_id.x"]
# 
# ggggg=subset(ICGC_PACA_AU_metaMatrix_duct_surv, icgc_donor_id.x %in% ICGC_PACA_AU_metaMatrix_duct_surv[duplicated(ICGC_PACA_AU_metaMatrix_duct_surv$icgc_donor_id.x),"icgc_donor_id.x"])

