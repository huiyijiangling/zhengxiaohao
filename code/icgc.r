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
# ICGC_PACA_CA = read_tsv(file="C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/ICGC/guanwang/PACA-CA/exp_seq.tsv")
# View(ICGC_PACA_CA)ICGC_PACA_CA
# system("cat C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/ICGC/guanwang/PACA-CA/exp_seq.tsv.gz |cut -f 1,3,8,9 >exp_seq_PRAD-CA_simplify.txt")
library(readr)
library(dplyr)
library(tidyverse)

donor= read_tsv(file="./ICGC/guanwang/PACA-CA/donor.tsv.gz",guess_max=min(1000000, Inf))#read_delim
specimen= read_tsv(file="./ICGC/guanwang/PACA-CA/specimen.tsv.gz",guess_max=min(1000000, Inf))#read_delim
sample=read_tsv(file="./ICGC/guanwang/PACA-CA/sample.tsv.gz",guess_max=min(1000000, Inf))#read_delim
therapy=read_tsv(file="./ICGC/guanwang/PACA-CA/donor_therapy.tsv.gz",guess_max=min(1000000, Inf))#read_delim
#TNM不一致按照高的来
parsing_failures <- problems(specimen)
parsing_failures
# kkkkkk<- list(donor, sample, specimen) %>% reduce(full_join, by = "icgc_donor_id") #这样其实产生了很多错误
# kkkkkk<- list(sample, specimen) %>% reduce(full_join, by = "icgc_specimen_id")
sample_specimen<-merge(sample,specimen,by="icgc_specimen_id",all=T)
sample_specimen_donor<-merge(sample_specimen,donor,by.x="icgc_donor_id.x",by.y="icgc_donor_id",all=T)
sample_specimen_donor=dplyr::rename(sample_specimen_donor,project_code.z=project_code.x,project_code.zz=project_code.y,submitted_donor_id.z=submitted_donor_id.x,submitted_donor_id.zz=submitted_donor_id.y)
sample_specimen_donor_therapy<-merge(sample_specimen_donor,therapy,by.x="icgc_donor_id.x",by.y="icgc_donor_id",all=T)
ICGC_PACA_CA_metaMatrix=sample_specimen_donor_therapy
table(duplicated(sample_specimen_donor_therapy$icgc_specimen_id))
table(duplicated(sample_specimen_donor_therapy$icgc_donor_id.x))
table(duplicated(sample_specimen_donor_therapy$icgc_sample_id))
# kkkkkk2<- list(sample, donor, specimen) %>% reduce(left_join, by = "icgc_donor_id")

#处理seq
ICGC_PACA_CA_seq = read_tsv(file="./ICGC/guanwang/PACA-CA/exp_seq.tsv.gz",guess_max=min(1000000, Inf))#read_delim
## 查看sample 数目
length(unique(ICGC_PACA_CA_seq$icgc_donor_id))    # [1] 91
length(unique(ICGC_PACA_CA_seq$icgc_specimen_id)) # [1] 92
length(unique(ICGC_PACA_CA_seq$icgc_sample_id)) # [1] 92 就用他
ICGC_PACA_CA_seq_only=unite(data=ICGC_PACA_CA_seq,col = "e",c("gene_id","icgc_sample_id"),remove = F)#"submitted_sample_id",
table(duplicated(ICGC_PACA_CA_seq_only$e))#有生物学重复
table(is.na(ICGC_PACA_CA_seq_only$normalized_read_count))#全都有
table(is.na(ICGC_PACA_CA_seq_only$raw_read_count))#全都有
ICGC_PACA_CA_seq_only=distinct(ICGC_PACA_CA_seq_only)
#这代码这么写的原因是删除重复的条目，只保留高表达量的条目
ICGC_PACA_CA_seq_only=ICGC_PACA_CA_seq_only[order(ICGC_PACA_CA_seq_only$e,ICGC_PACA_CA_seq_only$raw_read_count,ICGC_PACA_CA_seq_only$normalized_read_count,decreasing = T),]
ICGC_PACA_CA_seq_only=ICGC_PACA_CA_seq_only[!duplicated(ICGC_PACA_CA_seq_only$e),]
# ICGC_PACA_CA_seq_only_rep=ICGC_PACA_CA_seq_only[duplicated(ICGC_PACA_CA_seq_only$e),]#有重复的行，不限次数
# ICGC_PACA_CA_seq_only_rep=distinct(ICGC_PACA_CA_seq_only_rep,e)
# ICGC_PACA_CA_seq_only_rep=ICGC_PACA_CA_seq_only[ICGC_PACA_CA_seq_only$e %in% ICGC_PACA_CA_seq_only_rep$e,]
# ICGC_PACA_CA_seq_only_rep=ICGC_PACA_CA_seq_only_rep[order(ICGC_PACA_CA_seq_only_rep$e,ICGC_PACA_CA_seq_only_rep$raw_read_count,ICGC_PACA_CA_seq_only_rep$normalized_read_count,decreasing = T),]
#as.data.frame(table(ICGC_PACA_CA_seq_only_rep$e)>2)
#table(table(ICGC_PACA_CA_seq_only_rep$e)>2)#
table(table(ICGC_PACA_CA_seq_only$e))
table(duplicated(ICGC_PACA_CA_seq_only$e))
table(table(ICGC_PACA_CA_seq_only$icgc_donor_id))
table(ICGC_PACA_CA_seq_only$icgc_donor_id)
#icgc_sample_id,icgc_specimen_id是唯一的,但有时候如果生物重复过多，有可能出问题
#icgc_donor_id是重复的
#spread只能留一个列，多了否则就会报错
# gather 和spread 就是结构和合并的方法
ICGC_PACA_CA_seq_only_norm=ICGC_PACA_CA_seq_only[,c("normalized_read_count","gene_id","icgc_sample_id")]
ICGC_PACA_CA_seq_only_raw=ICGC_PACA_CA_seq_only[,c("raw_read_count","gene_id","icgc_sample_id")]
#
ICGC_PACA_CA_seq_only_norm <- ICGC_PACA_CA_seq_only_norm %>% group_by(gene_id) %>% spread(icgc_sample_id,normalized_read_count)
ICGC_PACA_CA_seq_only_norm=as.data.frame(ICGC_PACA_CA_seq_only_norm)
rownames(ICGC_PACA_CA_seq_only_norm)=ICGC_PACA_CA_seq_only_norm$gene_id
ICGC_PACA_CA_seq_only_norm=select(ICGC_PACA_CA_seq_only_norm,-"gene_id")
ICGC_PACA_CA_seq_only_raw <- ICGC_PACA_CA_seq_only_raw %>% group_by(gene_id)%>% spread(icgc_sample_id,raw_read_count)
ICGC_PACA_CA_seq_only_raw=as.data.frame(ICGC_PACA_CA_seq_only_raw)
rownames(ICGC_PACA_CA_seq_only_raw)=ICGC_PACA_CA_seq_only_raw$gene_id
ICGC_PACA_CA_seq_only_raw=select(ICGC_PACA_CA_seq_only_raw,-"gene_id")
if(F){
  #处理array
  ICGC_PACA_CA_array = read_tsv(file="./ICGC/guanwang/PACA-CA/exp_array.PACA-CA.tsv.gz")#read_delim
  ## 查看sample 数目
  length(unique(ICGC_PACA_CA_array$icgc_donor_id))    # [1] 91
  length(unique(ICGC_PACA_CA_array$icgc_specimen_id)) # [1] 92
  length(unique(ICGC_PACA_CA_array$icgc_sample_id)) # [1] 92 就用他
  #icgc_sample_id,icgc_specimen_id是唯一的
  #icgc_donor_id是重复的
  #spread只能留一个列，多了否则就会报错
  # gather 和spread 就是结构和合并的方法
  ICGC_PACA_CA_array_norm=ICGC_PACA_CA_array[,c("normalized_expression_value","gene_id","icgc_sample_id")]
  #
  ICGC_PACA_CA_array_norm <- ICGC_PACA_CA_array_norm %>% group_by(gene_id)  %>% spread(icgc_sample_id,normalized_expression_value)
  ICGC_PACA_CA_array_norm=as.data.frame(ICGC_PACA_CA_array_norm)
  rownames(ICGC_PACA_CA_array_norm)=ICGC_PACA_CA_array_norm$gene_id
  ICGC_PACA_CA_array_norm=select(ICGC_PACA_CA_array_norm,-"gene_id")#GSE36924 可以还原
}
write.csv(ICGC_PACA_CA_metaMatrix,file="ICGC_PACA_CA_metaMatrix.csv",quote = T)
ICGC_PACA_CA_seq_norm=ICGC_PACA_CA_seq_only_norm
ICGC_PACA_CA_seq_raw=ICGC_PACA_CA_seq_only_raw
ICGC_PACA_CA_seq=ICGC_PACA_CA_seq_only
#####################################
ICGC_PACA_CA_seq_log_filter=ICGC_PACA_CA_seq_raw
ICGC_PACA_CA_seq_log_filter[is.na(ICGC_PACA_CA_seq_log_filter)] <- 0# 我观察的要差就差一行整整,直接干掉
library(edgeR)
#生成not filter矩阵，选择F 20211206生成可cluster的矩阵
if(T){
  exprSet=ICGC_PACA_CA_seq_log_filter

  exprSet_rownames <- rownames(exprSet)
  exprSet_colnames <- colnames(exprSet)
  exprSet_quant <- preprocessCore::normalize.quantiles(
    as.matrix(exprSet))
  rownames(exprSet_quant) <- exprSet_rownames
  colnames(exprSet_quant) <- exprSet_colnames
  ICGC_PACA_CA_seq_log_all=log2(exprSet_quant+1)
}
#继续生成filter的矩阵
if(T){
expr = DGEList(counts = ICGC_PACA_CA_seq_log_filter)
keepALL <- rowSums(cpm(expr) > 1) >= 0.5*ncol(ICGC_PACA_CA_seq_log_filter)
nGenes <- as.numeric(summary(keepALL)[2]) + 
  as.numeric(summary(keepALL)[3])
nKeep <- summary(keepALL)[3]
cat (paste('Total Number of genes: ', nGenes, '\n', sep=''))
cat (paste('Number of genes for downstream analysis: ', nKeep, 
           '\n', sep=''))
exprSet=ICGC_PACA_CA_seq_log_filter[keepALL,]

# exprSet=rnaCounts_Pancreas[,group_list=="tumor"][keepALL,]
# dds <- dds[keep, ]
# dge <- dge[keep,,keep.lib.sizes = TRUE]
exprSet_rownames <- rownames(exprSet)
exprSet_colnames <- colnames(exprSet)
exprSet_quant <- preprocessCore::normalize.quantiles(
  as.matrix(exprSet))
rownames(exprSet_quant) <- exprSet_rownames
colnames(exprSet_quant) <- exprSet_colnames
ICGC_PACA_CA_seq_log_filter=log2(exprSet_quant+1)
}
save(ICGC_PACA_CA_seq_log_all,ICGC_PACA_CA_seq_log_filter,ICGC_PACA_CA_seq_norm,ICGC_PACA_CA_seq_raw,ICGC_PACA_CA_seq,ICGC_PACA_CA_metaMatrix,file = "ICGC_PACA_CA_expr.Rdata")
load("ICGC_PACA_CA_expr.Rdata")
#去掉杂质
#去掉经过术前治疗病历 保留 无治疗 手术 na 
table(ICGC_PACA_CA_metaMatrix$specimen_donor_treatment_type,useNA = "ifany")
table(ICGC_PACA_CA_metaMatrix$specimen_donor_treatment_type_other,useNA = "ifany")
ICGC_PACA_CA_metaMatrix=subset(ICGC_PACA_CA_metaMatrix,specimen_donor_treatment_type %in% c(NA,"surgery","no treatment"))
#duct
table(ICGC_PACA_CA_metaMatrix$specimen_type,useNA = "ifany")
table(ICGC_PACA_CA_metaMatrix$specimen_type_other,useNA = "ifany")
table(ICGC_PACA_CA_metaMatrix$tumour_histological_type,useNA = "ifany")
ICGC_PACA_CA_metaMatrix_ductNA=subset(ICGC_PACA_CA_metaMatrix,tumour_histological_type %in% c(NA,"8500/3"))
ICGC_PACA_CA_metaMatrix_ductNA_TAN=subset(ICGC_PACA_CA_metaMatrix_ductNA,specimen_type %in% c("Normal - solid tissue","Normal - tissue adjacent to primary","Primary tumour - solid tissue","Primary tumour - other"))#不能房间来NA，不是就不是

ICGC_PACA_CA_metaMatrix_ductNA_T=subset(ICGC_PACA_CA_metaMatrix_ductNA_TAN,specimen_type %in% c("Primary tumour - solid tissue","Primary tumour - other"))
ICGC_PACA_CA_metaMatrix_ductNA_A=subset(ICGC_PACA_CA_metaMatrix_ductNA_TAN,specimen_type %in% c("Normal - tissue adjacent to primary"))
ICGC_PACA_CA_metaMatrix_ductNA_N=subset(ICGC_PACA_CA_metaMatrix_ductNA_TAN,specimen_type %in% c("Normal - solid tissue"))
table(ICGC_PACA_CA_metaMatrix_ductNA_T$tumour_histological_type,useNA = "ifany")
table(ICGC_PACA_CA_metaMatrix_ductNA_A$tumour_histological_type,useNA = "ifany")
table(ICGC_PACA_CA_metaMatrix_ductNA_N$tumour_histological_type,useNA = "ifany")
#
ICGC_PACA_CA_metaMatrix_duct_T=subset(ICGC_PACA_CA_metaMatrix_ductNA_T,tumour_histological_type %in% c("8500/3"))
table(ICGC_PACA_CA_metaMatrix_duct_T$tumour_histological_type,useNA = "ifany")
#
# table(ICGC_PACA_CA_metaMatrix_duct_T$donor_vital_status)                   # counts the "unusual" NA 等于  table(ICGC_PACA_CA_metaMatrix_duct_T$donor_vital_status, exclude = NA)     # counts none
# table(ICGC_PACA_CA_metaMatrix_duct_T$donor_vital_status,useNA = "ifany")  # counts all three 等于  table(ICGC_PACA_CA_metaMatrix_duct_T$donor_vital_status, exclude = NULL)   #  (ditto)

# OS

ICGC_PACA_CA_metaMatrix_duct_T$OS.time=apply(ICGC_PACA_CA_metaMatrix_duct_T[,c("donor_relapse_interval","donor_interval_of_last_followup","donor_survival_time")],1,max)
ICGC_PACA_CA_metaMatrix_duct_T$OSbi=ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$donor_vital_status),NA,
                                           ifelse(ICGC_PACA_CA_metaMatrix_duct_T$donor_vital_status=="alive","0",
                                                  ifelse(ICGC_PACA_CA_metaMatrix_duct_T$donor_vital_status=="deceased","1","WRONG")))
ICGC_PACA_CA_metaMatrix_duct_T$OS=ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$OSbi)&is.na(ICGC_PACA_CA_metaMatrix_duct_T$OS.time),NA,
                                         ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$OSbi)&(!is.na(ICGC_PACA_CA_metaMatrix_duct_T$OS.time)),"0",
                                                ifelse(!is.na(ICGC_PACA_CA_metaMatrix_duct_T$OSbi),ICGC_PACA_CA_metaMatrix_duct_T$OSbi,"WRONG")))
# 某些人登记错的，rfi比OS还大# Relapse 其实是dfs和rfs不能混用的，不知道为啥说这个可以
#注意xx$xx时用==，na是不被考虑的
table(ICGC_PACA_CA_metaMatrix_duct_T$donor_relapse_interval,useNA = "ifany")#"stable"|"complete remission"|"no evidence of disease"
table(ICGC_PACA_CA_metaMatrix_duct_T$disease_status_last_followup,useNA = "ifany")#complete remission no evidence of disease partial remission progression  relapse  NA
#icgc这里最大的问题是没有 r0的数据。按照官网的说法，rfi等于dfi，但是cell2018的dfi及其严格，我更加倾向于这里面是包含r1结果的
#dfs也接受r1
#icgc的rfi属于宽指标，可以转换为dfs
#df严格要求，随访必须在90天以上
#RFI不需要分RFIbi，我算都一样
ICGC_PACA_CA_metaMatrix_duct_T$RFI=ifelse(!is.na(ICGC_PACA_CA_metaMatrix_duct_T$donor_relapse_interval),"1",
                                          ifelse(!is.na(ICGC_PACA_CA_metaMatrix_duct_T$donor_relapse_type),"1",
                                                 ifelse(ICGC_PACA_CA_metaMatrix_duct_T$disease_status_last_followup %in% c("relapse",'progression','partial remission'),"1",
                                                        ifelse(ICGC_PACA_CA_metaMatrix_duct_T$disease_status_last_followup %in% c("complete remission",'no evidence of disease','stable'),"0",
                                                               ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$disease_status_last_followup),NA,"Wrong")))))
ICGC_PACA_CA_metaMatrix_duct_T$RFI.time=ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$RFI),NA,
                                               ifelse(ICGC_PACA_CA_metaMatrix_duct_T$RFI=="0",ICGC_PACA_CA_metaMatrix_duct_T$OS.time,
                                                      ifelse(ICGC_PACA_CA_metaMatrix_duct_T$RFI=="1",ICGC_PACA_CA_metaMatrix_duct_T$donor_relapse_interval,"WRONG")))
ICGC_PACA_CA_metaMatrix_duct_T$DFS=ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$RFI),NA,
                                          ifelse(ICGC_PACA_CA_metaMatrix_duct_T$RFI=="0"&ICGC_PACA_CA_metaMatrix_duct_T$OS=="1","1",
                                                 ifelse(!is.na(ICGC_PACA_CA_metaMatrix_duct_T$RFI),ICGC_PACA_CA_metaMatrix_duct_T$RFI,"Wrong")))
ICGC_PACA_CA_metaMatrix_duct_T$DFS.time=ICGC_PACA_CA_metaMatrix_duct_T$RFI.time

table(ICGC_PACA_CA_metaMatrix_duct_T$DFS,useNA = "ifany")
table(ICGC_PACA_CA_metaMatrix_duct_T$DFS.time,useNA = "ifany")# Relapse type
table(ICGC_PACA_CA_metaMatrix_duct_T$donor_relapse_type,useNA = "ifany")# Relapse type
# DFI 去掉m1，去掉，没分期，去掉随访不足3月
# as.numeric(ICGC_PACA_CA_metaMatrix_duct_T$DFS.time))

#理论上应该是 给予一个确定值，如果这次是unknown就找上一个吗，不能直接只用NA，因为na得不出信息，截尾值也要知道从哪里截值才行。否则就是无意义信息。
#check stage info 最大作用就是直接干掉所有四期,保证术前的完整切除
#对于很多分析来说，4期是不纳入的
# relapse 取决于4期是否能完整切除 大多数不行
#DFI中是不纳入4期的。rfs也是。所有要求完整切除肿瘤的生存分析均有此要求
#请在这里修改，不要直接修改DFIOSPFS的代码
table(ICGC_PACA_CA_metaMatrix_duct_T$donor_tumour_stage_at_diagnosis,useNA = "ifany")# Relapse type
table(ICGC_PACA_CA_metaMatrix_duct_T$donor_tumour_stage_at_diagnosis_supplemental,useNA = "ifany")# Relapse type
table(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage,useNA = "ifany")# Relapse type
table(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage_supplemental,useNA = "ifany")# Relapse type
table(grepl("m1",ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage,ignore.case=T))
#具体修改
ICGC_PACA_CA_metaMatrix_duct_T_seestage=distinct(ICGC_PACA_CA_metaMatrix_duct_T[is.na(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage),c("icgc_sample_id","tumour_stage","tumour_stage_supplemental","donor_tumour_stage_at_diagnosis_supplemental","donor_tumour_stage_at_diagnosis")])
ICGC_PACA_CA_metaMatrix_duct_T$TNM=ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage
ICGC_PACA_CA_metaMatrix_duct_T[ICGC_PACA_CA_metaMatrix_duct_T$icgc_sample_id %in% "SA533649",]$TNM <- ICGC_PACA_CA_metaMatrix_duct_T[ICGC_PACA_CA_metaMatrix_duct_T$icgc_sample_id %in% "SA533649",]$donor_tumour_stage_at_diagnosis_supplemental
ICGC_PACA_CA_metaMatrix_duct_T$stage123=ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$donor_tumour_stage_at_diagnosis_supplemental)& is.na(ICGC_PACA_CA_metaMatrix_duct_T$donor_tumour_stage_at_diagnosis)&is.na(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage_supplemental)&is.na(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage),NA,
                                               ifelse(grepl("m1",ICGC_PACA_CA_metaMatrix_duct_T$TNM,ignore.case=T),"no",
                                                      ifelse(grepl("iv",ICGC_PACA_CA_metaMatrix_duct_T$TNM,ignore.case=T),"no","yes")))
if(F){
  #干掉m1
  ICGC_PACA_CA_metaMatrix_duct_T$DFS.time=ifelse(ICGC_PACA_CA_metaMatrix_duct_T$stage123=="no",NA,as.numeric(ICGC_PACA_CA_metaMatrix_duct_T$DFS.time))
  ICGC_PACA_CA_metaMatrix_duct_T$DFS=ifelse(ICGC_PACA_CA_metaMatrix_duct_T$stage123=="no",NA,ICGC_PACA_CA_metaMatrix_duct_T$DFS)
  
  ICGC_PACA_CA_metaMatrix_duct_T$RFI.time=ifelse(ICGC_PACA_CA_metaMatrix_duct_T$stage123=="no",NA,as.numeric(ICGC_PACA_CA_metaMatrix_duct_T$RFI.time))
  ICGC_PACA_CA_metaMatrix_duct_T$RFI=ifelse(ICGC_PACA_CA_metaMatrix_duct_T$stage123=="no",NA,ICGC_PACA_CA_metaMatrix_duct_T$RFI)
}
if(F){
  # 理论上还要干掉90天
ICGC_PACA_CA_metaMatrix_duct_T$DFI=ifelse(ICGC_PACA_CA_metaMatrix_duct_T$OS.time>=90,ICGC_PACA_CA_metaMatrix_duct_T$DFS,NA)
ICGC_PACA_CA_metaMatrix_duct_T$DFI.time=ifelse(ICGC_PACA_CA_metaMatrix_duct_T$OS.time>=90,ICGC_PACA_CA_metaMatrix_duct_T$DFS.time,NA)
#不建议直接套
}

# 理论上还要干掉r0
#不建议放这里处理
# ghjkl=subset(ICGC_PACA_CA_metaMatrix_duct_T,grepl("m1",ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage,ignore.case=T))
ICGC_PACA_CA_metaMatrix_duct_surv=ICGC_PACA_CA_metaMatrix_duct_T[,c("icgc_donor_id.x","DFS","DFS.time","OS","OS.time","RFI","RFI.time","stage123")]#"tumour_stage",
ICGC_PACA_CA_metaMatrix_duct_surv=dplyr::distinct(ICGC_PACA_CA_metaMatrix_duct_surv)
# tapply(ppp$OS.time,ppp$specimen_donor_treatment_type,median)
# qqq=ICGC_PACA_CA_metaMatrix_duct_T[grep("mx",ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage,ignore.case=T),]

# ICGC_PACA_CA_seq_norm,ICGCAU_PACA_CA_seq_raw,
library(AnnoProbe)
probe2gene=unique(ICGC_PACA_CA_metaMatrix_duct_T[,c("icgc_sample_id","icgc_donor_id.x")])#soft
colnames(probe2gene)=c("probe_id","symbol")

sampletodonor <- function(probes_expr){
  probes_expr=as.data.frame(t(probes_expr))
  genes_expr <- filterEM(probes_expr,probe2gene)
  genes_expr=as.data.frame(t(genes_expr))
  return(genes_expr)
}
ICGC_PACA_CA_seq_norm_dornor=sampletodonor(ICGC_PACA_CA_seq_norm)
ICGC_PACA_CA_seq_raw_dornor=sampletodonor(ICGC_PACA_CA_seq_raw)
ICGC_PACA_CA_seq_log_filter_dornor=sampletodonor(ICGC_PACA_CA_seq_log_filter)
ICGC_PACA_CA_seq_log_all_dornor=sampletodonor(ICGC_PACA_CA_seq_log_all)
#95
save(ICGC_PACA_CA_seq_log_all_dornor,ICGC_PACA_CA_seq_log_filter_dornor,ICGC_PACA_CA_seq_norm_dornor,ICGC_PACA_CA_seq_raw_dornor,ICGC_PACA_CA_metaMatrix_duct_surv,ICGC_PACA_CA_metaMatrix_duct_T,file = "ICGC_PACA_CA_metaMatrix_duct_surv.Rdata")
load("ICGC_PACA_CA_metaMatrix_duct_surv.Rdata")
write.csv(ICGC_PACA_CA_metaMatrix_duct_surv,file="ICGC_PACA_CA_metaMatrix_duct_surv.csv",quote = T)
write.csv(ICGC_PACA_CA_metaMatrix_duct_T,file="ICGC_PACA_CA_metaMatrix_T.csv",quote = T)
