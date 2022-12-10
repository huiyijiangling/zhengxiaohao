#第一步
#我认为实际上应该从tcga+geo validation用icgc
#先内部z-socre 在用sva
# v2 v7 交换

library(lars) 
library(glmnet) 
library(survival)
library(glmnet)
library(survival)
library(survminer)

  #注意计算时我直接使用OS置换了dfi
  library(readxl)
  # phe_gdc=read.csv(file="PAAD_ductal_clinical_gdc.csv“)
  #1
  phe_cellPAAD=read_excel("cellPAAD.xlsx",sheet = 1)
  phe_cellPAAD_OS01=phe_cellPAAD[phe_cellPAAD$OS %in% c(0,1),]
  dim(phe_cellPAAD_OS01)
  # phe_cellPAAD_OS01=as.matrix(phe_cellPAAD_OS01)
  # phe_cellPAAD_OS01[which(phe_cellPAAD_OS01=="NA")]<-NA
  # phe_cellPAAD_OS01=as.data.frame(phe_cellPAAD_OS01)
  
  if(T){
    #2
    phe_tcga_new=read_excel("PAAD_clinical_gdc.xlsx",sheet = 1,na="NA")#这个数据集直接删除好了。只有174个胰腺癌pancreatic cancer 不含有神经内，但有4个不正常表达
    phe_tcga_new=as.data.frame(phe_tcga_new)
    rownames(phe_tcga_new)=phe_tcga_new[,1]
    phe_tcga_new=phe_tcga_new[,-1]
    # phe_tcga_new[phe_tcga_new=="NA"]
    if(F){
      # phe_tcga_new=as.matrix(phe_tcga_new)
      # phe_tcga_new[which(phe_tcga_new=="NA")]<-NA
      # phe_tcga_new=as.data.frame(phe_tcga_new)
    }
    phe_tcga_new$anatomic_neoplasm_subdivision=ifelse(phe_tcga_new$anatomic_neoplasm_subdivision=="Other (please specify)",phe_tcga_new$anatomic_neoplasm_subdivision_other,phe_tcga_new$anatomic_neoplasm_subdivision)
    phe_tcga_new$histological_type=ifelse(phe_tcga_new$histological_type=="Pancreas-Adenocarcinoma-Other Subtype",phe_tcga_new$histological_type_other,phe_tcga_new$histological_type)
    # phe_tcga_new=phe_tcga_new[,c("bcr_patient_barcode","residual_tumor")]
    phe_cellPAAD_OS01=merge(phe_cellPAAD_OS01,phe_tcga_new,by = "bcr_patient_barcode")#TCGA-RL-AAAS 随访时间太短
    #who 认为谁是胰腺癌 2019
    #who 认为谁是胰腺癌 2010
    # 我的想法只删除 endocrine
    unique(phe_cellPAAD_OS01$histological_type.x)
    # [1] "Pancreas-Adenocarcinoma-Other Subtype"            "Pancreas-Adenocarcinoma Ductal Type"              "Pancreas-Undifferentiated Carcinoma"             
    # [4] "Pancreas-Colloid (mucinous non-cystic) Carcinoma" "NA"     
    unique(phe_cellPAAD_OS01$histological_type.y)
    # [1] "Pancreas-Adenocarcinoma-Other Subtype"            "Pancreas-Adenocarcinoma Ductal Type"              "Pancreas-Undifferentiated Carcinoma"             
    # [4] "Pancreas-Colloid (mucinous non-cystic) Carcinoma" "[Discrepancy]"    
    unique(phe_cellPAAD_OS01$histological_type_other)
    # [1] "invasive adenocarcinoma"                                                  "invasive, well-differentiated"                                           
    # [3] "NA"                                                                       "Invasive adenocarcinoma"                                                 
    # [5] "Poorly differentiated adenocarcinoma"                                     "Neuroendocrine"                                                          
    # [7] "NEUROENDOCRINE CARCINOMA NOS"                                             "82463 NEUROENDOCRINE CARCINOMA NOS"                                      
    # [9] "NEUROENDOCRINE CARCINOMA"                                                 "Adenocarcinoma, NOS"                                                     
    # [11] "Poorly differentiated pancreatic adenocarcinoma"                          "not specified"                                                           
    # [13] "Intraductal tubulopapillary neoplasm"                                     "ductal and micropapillary"                                               
    # [15] "Adenocarcinoma- NOS"                                                      "moderately differentiated ductal adenocarcinoma 60% + neuroendocrine 40%"
    table(phe_cellPAAD_OS01$histological_type.x=="[Discrepancy]")#	TCGA-YH-A8SY-01 这个人不详 删除
    phe_cellPAAD_OS01=subset(phe_cellPAAD_OS01,histological_type.x!="[Discrepancy]")
    table(grepl("NEUROENDOCRINE",phe_cellPAAD_OS01$histological_type_other,ignore.case=T))
    phe_cellPAAD_OS01=subset(phe_cellPAAD_OS01,!grepl("NEUROENDOCRINE",histological_type_other,ignore.case=T))
    #其实na和not specific 也要删
    # phe_cellPAAD_OS01=phe_cellPAAD_OS01[phe_cellPAAD_OS01$residual_tumor.y %in% c("R0","R1"),]
    # phe_cellPAAD_OS01$stage123=ifelse(grepl("m1",phe_cellPAAD_OS01$stage_event,ignore.case=T),"no",
    # ifelse(grepl("iv",phe_cellPAAD_OS01$pathologic_stage,ignore.case=T),"no","yes"))
    # phe_cellPAAD_OS01=subset(phe_cellPAAD_OS01,stage123=="yes")
  }
  
  if(F){
    library(tableone)
    
    phe_tcga_new$age65=ifelse(phe_tcga_new$age_at_initial_pathologic_diagnosis>=65,1,0)
    phe_tcga_new$nodepositive=ifelse(phe_tcga_new$number_of_lymphnodes_positive_by_he>0,1,0)
    phe_tcga_new$`Pathological stage (AJCC)`=ifelse(phe_tcga_new$pathologic_stage=="Stage IV","IV",
                                                    ifelse(is.na(phe_tcga_new$pathologic_stage),NA,
                                                           ifelse(!is.na(phe_tcga_new$pathologic_stage),"I-III",
                                                                  "Wrong")
                                                    )
    )
    phe_tcga_new$Differentiation=ifelse(phe_tcga_new$neoplasm_histologic_grade%in%c("G1","G2"),"Well or moderately differentiated",
                                        ifelse(phe_tcga_new$neoplasm_histologic_grade%in%c("G3","G4"),"Poorly differentiated",
                                               ifelse(phe_tcga_new$neoplasm_histologic_grade%in%c("GX"),NA,
                                                      ifelse(is.na(phe_tcga_new$neoplasm_histologic_grade),NA,"Wrong")
                                               )
                                        )
    )
    
    # ifelse(is.na(phe_tcga_new$residual_tumor),NA,
    phe_tcga_new$`Resection margin`=ifelse(phe_tcga_new$residual_tumor%in%c("R2","R1"),"Positive",
                                           ifelse(phe_tcga_new$residual_tumor%in%c("R0"),"Negative",
                                                  ifelse(phe_tcga_new$residual_tumor%in%c("RX"),NA,
                                                         ifelse(is.na(phe_tcga_new$residual_tumor),NA,"Wrong")
                                                  )
                                           )
    )
    phe_tcga_new$`Tumor location`=ifelse(phe_tcga_new$anatomic_neoplasm_subdivision %in% c("Pancreas Head and Body","Duct","Uncinate process","Head & Body of Pancreas","Head of Pancreas","Head and body","Head & body of pancreas","Head & body of pancreas","Head & Body of Pancreas"),"Head $ Neck",ifelse(phe_tcga_new$anatomic_neoplasm_subdivision %in% c("distal pancreas","Body of Pancreas","Body and Tail","Body & Tail of Pancreas","Tail of Pancreas"),"Body & Tail",ifelse(phe_tcga_new$anatomic_neoplasm_subdivision %in% c("overlapping parts of pancreas"),"Multiple",ifelse(is.na(phe_tcga_new$anatomic_neoplasm_subdivision),NA,"Wrong"))))
    #不需要先干掉na
    table(phe_tcga_new$neoplasm_histologic_grade,useNA = "ifany")
    # varsToFactor <- c("status","trt","ascites","hepato","spiders","edema","stage")
    # pbc[varsToFactor] <- lapply(pbc[varsToFactor], factor)
    ## Create Table 1 stratified by trt
    tableOne <- CreateTableOne(vars =  c("age65","gender","Tumor location","nodepositive","Pathological stage (AJCC)","Resection margin","Differentiation"),factorVars=c("age65","nodepositive"), data = phe_tcga_new, includeNA = F, addOverall = T)
    summary(tableOne)
    tableOne1=print(tableOne,formatOptions = list(big.mark = ","), quote = F, noSpaces = T, printToggle = T,missing=T,explain=F)
    write.csv(tableOne1,quote=T,file="tableOne1.csv")
    tableOne2=print(tableOne,  formatOptions = list(big.mark = ","), showAllLevels = TRUE, quote = F, noSpaces = T,missing=T, printToggle = T,explain=F)
    write.csv(tableOne2,quote=T,file="tableOne2.csv")
    ## Just typing the object name will invoke the print.TableOne method
    
  }
  #注意一定是用rnaExpr(已经配平的)，而不是用raw_count,raw count 没有校正。
  # metaMatrix.RNA=subset(metaMatrix.RNA,!sample %in%c("TCGA-3A-A9IV-01","TCGA-3A-A9IN-01","TCGA-3A-A9IL-01","TCGA-3A-A9IL-01","TCGA-2L-AAQM-01","TCGA-3A-A9IO-01","TCGA-3A-A9IR-01","TCGA-3A-A9IJ-01","TCGA-3A-A9IS-01","TCGA-YH-A8SY-01","TCGA-H8-A6C1-01"))#
  #其实还有
  #TCGA-HZ-A9TJ-06A	是一个meta标本已经删除
  #"TCGA-F2-6880-01"（表达谱mirrna不在pca内）
  #,"TCGA-H6-A45N-11","TCGA-2J-AABP-01" 表达谱rna有问题
  #TCGA-YH-A8SY-01 是NA
  phe=metaMatrix.RNA[substr(rownames(metaMatrix.RNA),14,15)=="01",]
  phe$event=ifelse(phe$vital_status=='Alive',0,1)
  phe$time=as.numeric(ifelse(phe$vital_status=='Alive',phe$days_to_last_follow_up,phe$days_to_death))
  #4 删除随访NA 
  dim(phe)
  # phe=as.matrix(phe)
  # phe[which(phe=="NA")]<-NA
  # phe=as.data.frame(phe)
  phe=subset(phe,!is.na(phe$time))
  #4 删除随访太少30，但是另一个研究直接写的1，怎么办滤过不掉啊
  dim(phe)
  # phe$time=phe$time/30
  phe$time=as.numeric(phe$time)
  phe=phe[phe$time>=30,]
  dim(phe)
  # phe=phe[-which(phe$time<36 & phe$event=="0"),]#不要随便删病例，对制造模型没有好处
  phe=merge(phe,phe_cellPAAD_OS01,by.x = "patient",by.y = "bcr_patient_barcode")#TCGA-RL-AAAS 随访时间太短
  is.na(phe$PFI)
  phe$DFS=ifelse(phe$OS==1,1,phe$PFI)
  is.na(phe$DFS)
  phe$DFS.time=phe$PFI.time
  # write.csv(phe,file="phe.csv",quote=T)
  
  
  
  phe$sex=ifelse(is.na(phe$gender),NA,
                 ifelse(phe$gender=="female",2,
                        ifelse(phe$gender=="male",1,"Wrong")))
  phe$fhpT=ifelse(is.na(phe$pathologic_T),NA,
                  ifelse(stringr::str_split(phe$pathologic_T,"T",simplify = T)[,2]%in%c("1","2","3","4"),stringr::str_split(phe$pathologic_T,"T",simplify = T)[,2],"Wrong"))
  phe$fhpT01=ifelse(is.na(phe$fhpT),NA,
                    ifelse(phe$fhpT%in%c(1,2),0,
                           ifelse(phe$fhpT%in%c(3,4),1,"Wrong")))
  phe$dcN01=ifelse(is.na(phe$pathologic_N),NA,
                   ifelse(stringr::str_split(phe$pathologic_N,'N',simplify = T)[,2]%in%c("1","1a","1b"),1,
                          ifelse(stringr::str_split(phe$pathologic_N,'N',simplify = T)[,2]%in%c("0","X"),0,"Wrong")))
  
  phe$diffg=phe$histological_grade
  phe$diff01=ifelse(is.na(phe$histological_grade),NA,
                    ifelse(phe$histological_grade %in% c("GX","Gx"),NA,
                           ifelse(phe$histological_grade%in%c("G1","G2"),0,
                                  ifelse(phe$histological_grade%in%c("G3","G4"),1,phe$histological_grade))))
  phe$stage02a=ifelse(is.na(phe$pathologic_stage),NA,
                      ifelse(phe$pathologic_stage%in%c("Stage IV","Stage III","Stage IIB"),1,#4
                             ifelse(phe$pathologic_stage%in%c("Stage IA","Stage IB","Stage IIA"),0,"Wrong"
                             )
                      )
  )
  phe$margin0=ifelse(is.na(phe$residual_tumor.y)|phe$residual_tumor.y %in% c("Rx","RX"),NA,
                     ifelse(phe$residual_tumor.y %in% c("R1","R2"),1,
                            ifelse(phe$residual_tumor.y %in% c("R0"),0,"Wrong")))
  phe$margin01=ifelse(is.na(phe$residual_tumor.y)|phe$residual_tumor.y %in% c("Rx","RX"),NA,
                      ifelse(phe$residual_tumor.y %in% c("R2"),1,
                             ifelse(phe$residual_tumor.y %in% c("R0","R1"),0,"Wrong")))
  phe=subset(phe,select=c("sample","sex","fhpT","fhpT01","dcN01","diff01","stage02a","margin0","margin01","OS.time","OS","DFS.time","DFS"))
  phe=unique(phe)
  phe[phe$sample %in% phe[duplicated(phe$sample),]$sample,]
  
  
  
  
  
  library(edgeR)
  expr = DGEList(counts = rnaCounts)
  keepALL <- rowSums(cpm(expr) > 1) >= 0.5*ncol(rnaCounts)
  nGenes <- as.numeric(summary(keepALL)[2]) + 
    as.numeric(summary(keepALL)[3])
  nKeep <- summary(keepALL)[3]
  cat (paste('Total Number of genes: ', nGenes, '\n', sep=''))
  cat (paste('Number of genes for downstream analysis: ', nKeep, 
             '\n', sep=''))
  exprSet=rnaCounts[,substr(colnames(rnaCounts),14,15)=="01"][keepALL,]
  # exprSet=rnaCounts_Pancreas[,group_list=="tumor"][keepALL,]
  # dds <- dds[keep, ]
  # dge <- dge[keep,,keep.lib.sizes = TRUE]
  exprSet_rownames <- rownames(exprSet)
  exprSet_colnames <- colnames(exprSet)
  exprSet_quant <- preprocessCore::normalize.quantiles(
    as.matrix(exprSet))
  rownames(exprSet_quant) <- exprSet_rownames
  colnames(exprSet_quant) <- exprSet_colnames
  exprSet_quant=log2(exprSet_quant+1)
  #是不是应该 keepsize
  # if(T){
  #   load("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/ans.Rdata")
  #   library(stringr)
  #   lnchighall=unique(str_split(ansT3,"_",simplify = T)[,1])#ansALL
  #   mhighall=unique(str_split(ansT3,"_",simplify = T)[,2])#ansALL
  #   # ansT
  #   # ansN
  # }
  # #这里面建模时候就要考虑到overlap
  # 
  mhighall=rownames(exprSet_quant_with_clinical)
  #############################################
  if(T){
    load("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/EMTAB6134_after_bioc.Rdata")
    load("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/ICGC_PACA_AU_metaMatrix_duct_surv.Rdata")
    load("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/ICGC_PACA_CA_metaMatrix_duct_surv.Rdata")
    load("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/GSE71729_after_soft.Rdata")
    load("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/GSE21501_after_soft.Rdata")#external validation
    lnchighall_over=Reduce(intersect,list(lnchighall,rownames(ICGC_PACA_CA_seq_log_filter_dornor),rownames(ICGC_PACA_AU_array_norm_dornor),rownames(genes_expr_mean_EMTAB6134),rownames(genes_expr_mean_GSE71729_rna),rownames(genes_expr_mean_GSE21501_rna)))#,rownames(ICGC_PACA_AU_seq_log_filter_dornor)
    mhighall_over=Reduce(intersect,list(mhighall,rownames(ICGC_PACA_CA_seq_log_filter_dornor),rownames(ICGC_PACA_AU_array_norm_dornor),rownames(genes_expr_mean_EMTAB6134),rownames(genes_expr_mean_GSE71729_rna),rownames(genes_expr_mean_GSE21501_rna)))#,rownames(ICGC_PACA_AU_seq_log_filter_dornor)
  }
  
  #中心化前去掉不分析的样本
  phe=subset(phe,!sample %in%c("TCGA-F2-6880-01","TCGA-H6-A45N-11","TCGA-2J-AABP-01"))#pca 分析
  
  exprSet_quant_with_clinical=exprSet_quant[,which(colnames(exprSet_quant) %in% phe$sample)]#[mhighall_over,]#rownames(dePC)
  phe=phe[which(phe$sample %in% colnames(exprSet_quant_with_clinical)),]
  phe=phe[order(phe$sample,decreasing = F),]
  exprSet_quant_with_clinical=exprSet_quant_with_clinical[,order(colnames(exprSet_quant_with_clinical),decreasing = F)]

  
  
  library(lars) 
  library(glmnet) 
  library(survival)
  library(glmnet)
  library(survival)
  library(survminer)
  
  all(substring(colnames(exprSet_quant_with_clinical),1,12)==phe$patient)
  tcga_only_expr_mean=t(scale(t(exprSet_quant_with_clinical)))
  x=t(tcga_only_expr_mean)#as.matrix(t(exprSet_quant_with_clinical))
  # y <- Surv(phe$DFS.time, phe$DFS)#Surv(phe$time, phe$event)
  y <- Surv(phe$OS.time, phe$OS)#Surv(phe$time, phe$event)
  # fit the model
  x=t(tcga_only_expr_mean)
  dat=as.data.frame(x)
  dat$sample=rownames(dat)
  dat=merge(phe,dat,by="sample")
  library(survival)
  library(survminer)
  dat$OS=as.numeric(dat$OS)
  dat$OS.time=as.numeric(dat$OS.time)
  dat$DFS=as.numeric(dat$DFS)
  dat$DFS.time=as.numeric(dat$DFS.time)
  dat$event=as.numeric(dat$OS)
  dat$time=as.numeric(dat$OS.time)
  # mat.100 <- matrix(nrow= 100, ncol=2)
  # mat.out <- matrix(nrow= 100, ncol=7)
  }
  
  
  
  
  
  #F immune
  # immune_timer2=read.csv(file ="C:/Users/zxh/Desktop/R/ftoimmune/timer2/infiltration_estimation_for_tcga.csv.gz",header = T,stringsAsFactors = F )
  # 
  # rownames(immune_timer2)=immune_timer2$cell_type
  # immune_timer2$cell_type=NULL
  # immune_timer2=t(immune_timer2)
  # immune_timer2=as.data.frame(immune_timer2)
  # table(stringr::str_split(rownames(immune_timer2),"_",simplify = T)[,2])
  # # immune_timer2=subset(immune_timer2,stringr::str_split(rownames(immune_timer2),"_",simplify = T)[,2] %in%"XCELL")
  # # rownames(immune_timer2_CIBERSORT.ABS)=stringr::str_split(rownames(immune_timer2_CIBERSORT.ABS),"_",simplify = T)[,1]
  # # rownames(immune_timer2_CIBERSORT.ABS)=gsub("\\.", "_", rownames(immune_timer2_CIBERSORT.ABS)) 
  # 
  # # immune_timer2=subset(immune_timer2,stringr::str_split(rownames(immune_timer2),"_",simplify = T)[,2] %in%"CIBERSORT.ABS")
  # # rownames(immune_timer2)=stringr::str_split(rownames(immune_timer2),"_",simplify = T)[,1]
  # rownames(immune_timer2)=gsub("\\.", "_", rownames(immune_timer2))
  # 
  # immune_timer2=subset(immune_timer2,select = phe$sample)
  # write.csv(v0,"v0.csv",quote = T)
  # if(F){
  #   #F merrrr=subset(merrrr,!merrrr$seriesall=="v2")#  merrrr=subset(merrrr,!merrrr$seriesall=="v1")
  #   # ui=subset(ui,select = merrrr$sample)
  #   
  #   v0=as.data.frame(t(immune_timer2))
  #   v0$sample=rownames(v0)
  #   v0=merge(phe,v0,by="sample")
  #   v0$event=as.numeric(v0$OS)
  #   v0$time=as.numeric(v0$OS.time)
  #   
  #   x=t(immune_timer2)
  #   
  #   y <- Surv(v0$OS.time, v0$OS)
  # }
  
  if(T){
    #没有age
    v2=t(genes_expr_mean_EMTAB6134)
    v2=as.data.frame(v2)
    v2$geo_accession=rownames(v2)
    phenoDat_EMTAB6134=subset(phenoDat_EMTAB6134,select=c('Source Name','Characteristics[sex]','Characteristics[tumor grading]','Characteristics[TNM tumour grading]','Characteristics[resection margin]','Characteristics[os.delay]','Characteristics[os.event]','Characteristics[dfs.delay]','Characteristics[dfs.event]','OS.time','OS','DFS.time','DFS'))
    # 'Characteristics[hightumcellclassif]','Characteristics[wholetumclassif]','Characteristics[average vaf]','Characteristics[krasmut]','Characteristics[tp53mut]','Characteristics[cdkn2amut]',
    colnames(phenoDat_EMTAB6134)=c("Source Name","sex","diffg","TN","margin","orginal.os.time","orginal.os.event","orginal.dfs.time","orginal.dfs.event","OS.time","OS","DFS.time","DFS")#"hightumcellclassif","wholetumclassif","averagevaf","krasmut","tp53mut","cdkn2amut"
    phenoDat_EMTAB6134=do::Replace(data=phenoDat_EMTAB6134, from=c("not available"),to=NA)
    phenoDat_EMTAB6134$fhpT=substr(phenoDat_EMTAB6134$TN,2,2)
    phenoDat_EMTAB6134$fhpT01=ifelse(phenoDat_EMTAB6134$fhpT%in%c(1,2),0,
                                     ifelse(phenoDat_EMTAB6134$fhpT%in%c(3,4),1,"Wrong"))
    phenoDat_EMTAB6134$dcN01=substr(phenoDat_EMTAB6134$TN,4,4)
    phenoDat_EMTAB6134$diffg=stringr::str_split(phenoDat_EMTAB6134$diffg," ",simplify = T)[,3]
    phenoDat_EMTAB6134$diff01=ifelse(phenoDat_EMTAB6134$diffg %in% c("G1","G2"),0,ifelse(phenoDat_EMTAB6134$diffg %in% c("G3","G4"),1,"Wrong"))
    phenoDat_EMTAB6134$sex=ifelse(is.na(phenoDat_EMTAB6134$sex),NA,
                                  ifelse(phenoDat_EMTAB6134$sex=="female",2,
                                         ifelse(phenoDat_EMTAB6134$sex=="male",1,"Wrong")))
    phenoDat_EMTAB6134$margin0=ifelse(is.na(phenoDat_EMTAB6134$margin),NA,
                                      ifelse(phenoDat_EMTAB6134$margin %in% c("resection margin R1","resection margin R2"),1,
                                             ifelse(phenoDat_EMTAB6134$margin %in% c("resection margin R0"),0,"Wrong")))
    phenoDat_EMTAB6134$margin01=ifelse(is.na(phenoDat_EMTAB6134$margin),NA,
                                       ifelse(phenoDat_EMTAB6134$margin %in% c("resection margin R1","resection margin R0"),0,
                                              ifelse(phenoDat_EMTAB6134$margin %in% c("resection margin R2"),1,"Wrong")))
    phenoDat_EMTAB6134$stage02a=ifelse(is.na(phenoDat_EMTAB6134$dcN01),NA,
                                       ifelse(phenoDat_EMTAB6134$dcN01==1|phenoDat_EMTAB6134$fhpT==4,1,
                                              ifelse(phenoDat_EMTAB6134$dcN01==0,0,"Wrong")))
    phenoDat_EMTAB6134=subset(phenoDat_EMTAB6134,select=c("Source Name","sex","fhpT","fhpT01","dcN01","diff01","stage02a","margin0","margin01","OS.time","OS","DFS.time","DFS"))
    v2=merge(phenoDat_EMTAB6134,v2,by.x="Source Name",by.y="geo_accession")
    v2$event=as.numeric(v2$OS)
    v2$OS.time=as.numeric(v2$OS.time)*30
    v2$OS.time=round(v2$OS.time,0)
    v2$time=round(v2$OS.time,0)
    ###################################################
    ICGC_PACA_CA_seq_log_filter_dornor_mean=t(scale(t(ICGC_PACA_CA_seq_log_filter_dornor)))
    v3=t(ICGC_PACA_CA_seq_log_filter_dornor_mean)
    v3=as.data.frame(v3)
    v3$icgc_donor_id.x=rownames(v3)
    dim(ICGC_PACA_CA_metaMatrix_duct_T)[1]==dim(ICGC_PACA_CA_metaMatrix_duct_surv)[1]#false
    ICGC_PACA_CA_metaMatrix_duct_T=subset(ICGC_PACA_CA_metaMatrix_duct_T,select=c("icgc_donor_id.x","tumour_grade","tumour_stage","study_donor_involved_in","donor_sex","donor_vital_status","disease_status_last_followup","donor_relapse_type","donor_age_at_diagnosis","donor_age_at_enrollment","donor_age_at_last_followup","donor_relapse_interval","donor_diagnosis_icd10","donor_tumour_stage_at_diagnosis","donor_survival_time","donor_interval_of_last_followup","prior_malignancy","cancer_type_prior_malignancy","cancer_history_first_degree_relative","project_code.y.1","submitted_donor_id.y.1","first_therapy_type","first_therapy_therapeutic_intent","first_therapy_start_interval","first_therapy_duration","first_therapy_response","second_therapy_type","second_therapy_therapeutic_intent","second_therapy_start_interval","second_therapy_duration","second_therapy_response","other_therapy","other_therapy_response","OS.time","OSbi","OS","RFI","RFI.time","DFS","DFS.time","TNM","stage123"))
    # write.csv(ICGC_PACA_CA_metaMatrix_duct_T,file = "ICGC_PACA_CA_metaMatrix_duct_T.csv",quote = T)
    ICGC_PACA_CA_metaMatrix_duct_T=unique(ICGC_PACA_CA_metaMatrix_duct_T)
    dim(ICGC_PACA_CA_metaMatrix_duct_T)[1]==dim(ICGC_PACA_CA_metaMatrix_duct_surv)[1]
    ICGC_PACA_CA_metaMatrix_duct_T$sex=ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$donor_sex),NA,
                                              ifelse(ICGC_PACA_CA_metaMatrix_duct_T$donor_sex=="female",2,
                                                     ifelse(ICGC_PACA_CA_metaMatrix_duct_T$donor_sex=="male",1,"Wrong")))
    ICGC_PACA_CA_metaMatrix_duct_T$fhpT=NA
    ICGC_PACA_CA_metaMatrix_duct_T$fhpT01=NA
    ICGC_PACA_CA_metaMatrix_duct_T$dcN01=NA
    ICGC_PACA_CA_metaMatrix_duct_T$diffg=ICGC_PACA_CA_metaMatrix_duct_T$tumour_grade
    ICGC_PACA_CA_metaMatrix_duct_T$diff01=ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$tumour_grade),NA,
                                                 ifelse(ICGC_PACA_CA_metaMatrix_duct_T$tumour_grade%in%c("Well differentiated","Moderately differentiated"),0,
                                                        ifelse(ICGC_PACA_CA_metaMatrix_duct_T$tumour_grade%in%c("Poorly differentiated","Undifferentiated","Moderate to Poor"),1,ICGC_PACA_CA_metaMatrix_duct_T$tumour_grade)))
    ICGC_PACA_CA_metaMatrix_duct_T$stage02a=ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage),NA,
                                                   ifelse(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage %in% c("IV","III","IIB"),1,
                                                          ifelse(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage %in% c("I","IA","IB","II","IIA"),0,"Wrong")))
    ICGC_PACA_CA_metaMatrix_duct_T$margin01=NA
    ICGC_PACA_CA_metaMatrix_duct_T$margin0=NA
    
    
    # ICGC_PACA_CA_metaMatrix_duct_T$stage24=ifelse(is.na(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage),NA,
    #                                                ifelse(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage %in% c("IV","III"),1,
    #                                                       ifelse(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage %in% c("I","IA","IB","II","IIA","IIB"),0,"Wrong")))
    # table(ICGC_PACA_CA_metaMatrix_duct_T$tumour_stage)
    ICGC_PACA_CA_metaMatrix_duct_T=subset(ICGC_PACA_CA_metaMatrix_duct_T,select=c("icgc_donor_id.x","sex","fhpT","fhpT01","dcN01","diff01","stage02a","margin0","margin01","OS.time","OS","DFS.time","DFS"))
    ICGC_PACA_CA_metaMatrix_duct_T=unique(ICGC_PACA_CA_metaMatrix_duct_T)
    v3=merge(ICGC_PACA_CA_metaMatrix_duct_T,v3,by="icgc_donor_id.x")#ICGC_PACA_CA_metaMatrix_duct_surv
    v3$event=as.numeric(v3$OS)
    v3$time=as.numeric(v3$OS.time)
    v3=v3[!is.na(v3$OS),]
    v3=v3[!is.na(v3$OS.time),]
    # v3$event=as.numeric(v3$DFS)
    # v3$time=as.numeric(v3$DFS.time)  
    v3[duplicated(v3$icgc_donor_id.x),"icgc_donor_id.x"]
    ############################################################
    ICGC_PACA_AU_array_norm_dornor_mean=t(scale(t(ICGC_PACA_AU_array_norm_dornor)))
    v4=t(ICGC_PACA_AU_array_norm_dornor_mean)
    v4=as.data.frame(v4)
    v4$icgc_donor_id.x=rownames(v4)
    
    ICGC_PACA_AU_metaMatrix_duct_T$sex=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$donor_sex),NA,
                                              ifelse(ICGC_PACA_AU_metaMatrix_duct_T$donor_sex=="female",2,
                                                     ifelse(ICGC_PACA_AU_metaMatrix_duct_T$donor_sex=="male",1,"Wrong")))
    ICGC_PACA_AU_metaMatrix_duct_T$fhpT=stringr::str_split(ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,'[TNM]',simplify = T)[,2]
    ICGC_PACA_AU_metaMatrix_duct_T$fhpT=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$fhpT)|ICGC_PACA_AU_metaMatrix_duct_T$fhpT==""|ICGC_PACA_AU_metaMatrix_duct_T$fhpT=="X"|ICGC_PACA_AU_metaMatrix_duct_T$fhpT==" ",NA,
                                               ifelse(ICGC_PACA_AU_metaMatrix_duct_T$fhpT%in%c(1,2,3,4),ICGC_PACA_AU_metaMatrix_duct_T$fhpT,"Wrong"))
    ICGC_PACA_AU_metaMatrix_duct_T$fhpT01=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$fhpT),NA,
                                                 ifelse(ICGC_PACA_AU_metaMatrix_duct_T$fhpT%in%c(1,2),0,
                                                        ifelse(ICGC_PACA_AU_metaMatrix_duct_T$fhpT%in%c(3,4),1,"Wrong")))
    ICGC_PACA_AU_metaMatrix_duct_T$dcN01=stringr::str_split(ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,'[TNM]',simplify = T)[,3]
    ICGC_PACA_AU_metaMatrix_duct_T$dcN01=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$dcN01)|ICGC_PACA_AU_metaMatrix_duct_T$dcN01==""|ICGC_PACA_AU_metaMatrix_duct_T$dcN01=="X"|ICGC_PACA_AU_metaMatrix_duct_T$dcN01==" ",NA,
                                                ifelse(ICGC_PACA_AU_metaMatrix_duct_T$dcN01%in%c("1","1a","1b"),1,
                                                       ifelse(ICGC_PACA_AU_metaMatrix_duct_T$dcN01==c("0"),0,"Wrong")))
    ICGC_PACA_AU_metaMatrix_duct_T$diffg=ICGC_PACA_AU_metaMatrix_duct_T$tumour_grade
    ICGC_PACA_AU_metaMatrix_duct_T$diff01=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$tumour_grade),NA,
                                                 ifelse(ICGC_PACA_AU_metaMatrix_duct_T$tumour_grade %in% c("X - Cannot be assessed"),NA,
                                                        ifelse(ICGC_PACA_AU_metaMatrix_duct_T$tumour_grade%in%c("1 - Well differentiated","2 - Moderately differentiated"),0,
                                                               ifelse(ICGC_PACA_AU_metaMatrix_duct_T$tumour_grade%in%c("3 - Poorly differentiated","4 - Undifferentiated"),1,ICGC_PACA_AU_metaMatrix_duct_T$tumour_grade))))
    # table(grepl("NX",ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,ignore.case=T))
    # table(grepl("TX",ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,ignore.case=T))
    # table(ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage[grepl("NX",ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,ignore.case=T)])
    # table(ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage[grepl("TX",ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,ignore.case=T)])
    ICGC_PACA_AU_metaMatrix_duct_T$stage02a=ifelse(is.na(ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage)|ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage=="TXNXMX",NA,
                                                   ifelse(grepl("M1",ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,ignore.case=T),1,#4
                                                          ifelse(grepl("T4",ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,ignore.case=T),1,#3
                                                                 ifelse(grepl("N1",ICGC_PACA_AU_metaMatrix_duct_T$tumour_stage,ignore.case=T),1,0)
                                                          )
                                                   )
    )
    ICGC_PACA_AU_metaMatrix_duct_T$margin0=NA
    ICGC_PACA_AU_metaMatrix_duct_T$margin01=NA
    ICGC_PACA_AU_metaMatrix_duct_T=subset(ICGC_PACA_AU_metaMatrix_duct_T,select=c("icgc_donor_id.x","sex","fhpT","fhpT01","dcN01","diff01","stage02a","margin0","margin01","OS.time","OS","DFS.time","DFS"))
    ICGC_PACA_AU_metaMatrix_duct_T=unique(ICGC_PACA_AU_metaMatrix_duct_T)
    ICGC_PACA_AU_metaMatrix_duct_T[ICGC_PACA_AU_metaMatrix_duct_T$icgc_donor_id.x %in% ICGC_PACA_AU_metaMatrix_duct_T[duplicated(ICGC_PACA_AU_metaMatrix_duct_T$icgc_donor_id.x),]$icgc_donor_id.x,]
    ICGC_PACA_AU_metaMatrix_duct_T=subset(ICGC_PACA_AU_metaMatrix_duct_T,
                                          !(ICGC_PACA_AU_metaMatrix_duct_T$icgc_donor_id.x=="DO32831"&ICGC_PACA_AU_metaMatrix_duct_T$dcN01=="0"))
    ICGC_PACA_AU_metaMatrix_duct_T=subset(ICGC_PACA_AU_metaMatrix_duct_T,
                                          !(ICGC_PACA_AU_metaMatrix_duct_T$icgc_donor_id.x=="DO32881"&ICGC_PACA_AU_metaMatrix_duct_T$diff01=="0")) 
    ICGC_PACA_AU_metaMatrix_duct_T=subset(ICGC_PACA_AU_metaMatrix_duct_T,
                                          !(ICGC_PACA_AU_metaMatrix_duct_T$icgc_donor_id.x=="DO33063"&ICGC_PACA_AU_metaMatrix_duct_T$fhpT=="2"))  
    v4=merge(ICGC_PACA_AU_metaMatrix_duct_T,v4,by="icgc_donor_id.x")#ICGC_PACA_AU_metaMatrix_duct_surv
    v4$event=as.numeric(v4$OS)
    v4$time=as.numeric(v4$OS.time)
    v4=v4[!is.na(v4$OS),]
    v4=v4[!is.na(v4$OS.time),]
    # v4$event=as.numeric(v4$DFS)
    # v4$time=as.numeric(v4$DFS.time)  
    v4[v4$icgc_donor_id.x%in%v4[duplicated(v4$icgc_donor_id.x),"icgc_donor_id.x"],][,1:12]
    if(F){
      ICGC_PACA_AU_seq_log_filter_dornor_mean=t(scale(t(ICGC_PACA_AU_seq_log_filter_dornor)))
      v5=t(ICGC_PACA_AU_seq_log_filter_dornor_mean)
      v5=as.data.frame(v5)
      v5$icgc_donor_id.x=rownames(v5)
      v5=merge(ICGC_PACA_AU_metaMatrix_duct_surv,v5,by="icgc_donor_id.x")
      v5$event=as.numeric(v5$OS)
      v5$time=as.numeric(v5$OS.time)
      v5=v5[!is.na(v5$OS),]
      v5=v5[!is.na(v5$OS.time),]
      # v5$event=as.numeric(v5$DFS)
      # v5$time=as.numeric(v5$DFS.time)
      v5[duplicated(v5$icgc_donor_id.x),"icgc_donor_id.x"]
      #v5 没有意义 因为大多数都没生存时间，而且都在array里有了
    }
    v6=t(genes_expr_mean_GSE71729_rna)
    v6=as.data.frame(v6)
    v6$geo_accession=rownames(v6)
    phenoDat_GSE71729_rna$sex=NA
    phenoDat_GSE71729_rna$fhpT=NA
    phenoDat_GSE71729_rna$fhpT01=NA
    phenoDat_GSE71729_rna$dcN01=NA
    phenoDat_GSE71729_rna$diffg=NA
    phenoDat_GSE71729_rna$diff01=NA
    phenoDat_GSE71729_rna$stage02a=NA
    phenoDat_GSE71729_rna$margin0=NA
    phenoDat_GSE71729_rna$margin01=NA
    phenoDat_GSE71729_rna$DFS=NA
    phenoDat_GSE71729_rna$DFS.time=NA  
    phenoDat_GSE71729_rna=subset(phenoDat_GSE71729_rna,select=c("geo_accession","sex","fhpT","fhpT01","dcN01","diff01","stage02a","margin0","margin01","OS.time","OS","DFS.time","DFS"))
    v6=merge(phenoDat_GSE71729_rna,v6,by.x="geo_accession",by.y="geo_accession")
    v6$event=as.numeric(v6$OS)
    v6$OS.time=as.numeric(v6$OS.time)*30
    v6$time=round(v6$OS.time,0)
    #这个连TN都没有我去
    ##########################################################
    v7=t(genes_expr_mean_GSE21501_rna)
    v7=as.data.frame(v7)
    v7$geo_accession=rownames(v7)
    phenoDat_GSE21501_rna$fhpT01=ifelse(is.na(phenoDat_GSE21501_rna$fhpT),NA,
                                        ifelse(phenoDat_GSE21501_rna$fhpT%in%c(1,2),0,
                                               ifelse(phenoDat_GSE21501_rna$fhpT%in%c(3,4),1,"Wrong")))
    phenoDat_GSE21501_rna=subset(phenoDat_GSE21501_rna,select=c("geo_accession","sex","fhpT","fhpT01","dcN01","diff01","stage02a","margin0","margin01","OS.time","OS","DFS.time","DFS"))
    v7=merge(phenoDat_GSE21501_rna,v7,by.x="geo_accession",by.y="geo_accession")
    v7$event=as.numeric(v7$OS)
    v7$OS.time=as.numeric(v7$OS.time)
    v7$time=round(v7$OS.time,0)
    #
    #理论上可以用sva校正的。
    #还有internal和external来说呢。
    #这个连TN都没有我去
    
  }
  # immune_timer2_CIBERSORT.ABS=subset(immune_timer2,stringr::str_split(rownames(immune_timer2),"_",simplify = T)[,2] %in%"CIBERSORT.ABS")
  # rownames(immune_timer2_CIBERSORT.ABS)=stringr::str_split(rownames(immune_timer2_CIBERSORT.ABS),"_",simplify = T)[,1]
  # rownames(immune_timer2_CIBERSORT.ABS)=gsub("\\.", "_", rownames(immune_timer2_CIBERSORT.ABS)) 
  # 
  # immune_timer2_CIBERSORT.ABS_dataset=subset(immune_timer2_CIBERSORT.ABS,select = phe$sample)
  
  ppp_FTO_2=read.csv("C:/Users/zxh/Desktop/R/ftoimmune/ppp_FTO_2.csv")
  load("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/ppp_icgc_ca_dornor_lose3.Rdata")
  
  if(T){

    v0=as.data.frame(t(TCGA_immune))
    v0$sample=rownames(v0)
    v0=merge(phe,v0,by="sample")
    v0$event=as.numeric(v0$OS)
    v0$time=as.numeric(v0$OS.time)
    x=t(subset(TCGA_immune,select=v0$sample,rownames(TCGA_immune)%in%ppp_FTO_2$X))
    
    y <- Surv(v0$OS.time, v0$OS)
  }
  
  colnames(x)=gsub(" ", "_", colnames(x))
  colnames(x)=gsub("\\+", "_", colnames(x))
  colnames(x)=gsub("-", "_", colnames(x))
  colnames(v0)=gsub(" ", "_", colnames(v0))
  colnames(v0)=gsub("\\+", "_", colnames(v0))
  colnames(v0)=gsub("-", "_", colnames(v0))
  rownames(icgc_CA_tpm_immune_dornor)=gsub(" ", "_", rownames(icgc_CA_tpm_immune_dornor))
  rownames(icgc_CA_tpm_immune_dornor)=gsub("\\+", "_", rownames(icgc_CA_tpm_immune_dornor))
  rownames(icgc_CA_tpm_immune_dornor)=gsub("-", "_", rownames(icgc_CA_tpm_immune_dornor))
  ppp_FTO_2$X=gsub(" ", "_", ppp_FTO_2$X)
  ppp_FTO_2$X=gsub("\\+", "_", ppp_FTO_2$X)
  ppp_FTO_2$X=gsub("-", "_", ppp_FTO_2$X)
  
  # ICGC_PACA_CA_metaMatrix_duct_T=subset(ICGC_PACA_CA_metaMatrix_duct_T,select=c("icgc_donor_id.x","sex","fhpT","fhpT01","dcN01","diff01","stage02a","margin0","margin01","OS.time","OS","DFS.time","DFS"))
  # ICGC_PACA_CA_metaMatrix_duct_T=unique(ICGC_PACA_CA_metaMatrix_duct_T)

  if(T){
    icgc_CA_tpm_immune_dornor=subset(icgc_CA_tpm_immune_dornor,rownames(icgc_CA_tpm_immune_dornor)%in%ppp_FTO_2$X)
    icgc_CA_tpm_immune_dornor=as.data.frame(t(icgc_CA_tpm_immune_dornor))
    icgc_CA_tpm_immune_dornor$icgc_donor_id.x=rownames(icgc_CA_tpm_immune_dornor)
    icgc_CA_tpm_immune_dornor=merge(ICGC_PACA_CA_metaMatrix_duct_T,icgc_CA_tpm_immune_dornor,by="icgc_donor_id.x",all = F)
    icgc_CA_tpm_immune_dornor$event=as.numeric(icgc_CA_tpm_immune_dornor$OS)
    icgc_CA_tpm_immune_dornor$time=as.numeric(icgc_CA_tpm_immune_dornor$OS.time)
    icgc_CA_tpm_immune_dornor=subset(icgc_CA_tpm_immune_dornor,!is.na(icgc_CA_tpm_immune_dornor$time))
  }


  if(T){
    FGT3 <- function(i){
      library(survival)
      set.seed(i)
      
      cv_fit <- glmnet::cv.glmnet(x, y, alpha=1,type.measure="C",standardize = FALSE,family = 'cox',nfolds =3)
      fit <- glmnet::glmnet(x=x, y=y, alpha = 1,family ='cox', lambda=cv_fit$lambda.1se)#lambda.1se
      choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
      if(length(choose_gene)!=0){s=as.formula(paste("Surv(time, event) ~" ,paste(choose_gene,collapse = "+")))
      model <- survival::coxph(formula=s, data =v0)
      #
      new_dat=v0
      new_dat$event=as.numeric(new_dat$event)
      new_dat$time=as.numeric(new_dat$time)
      fp <- predict(model,new_dat) ;
      library(Hmisc)
      options(scipen=200)
      cc0=1-with(new_dat,rcorr.cens(fp,Surv(time, event)))["C Index"]

      mat.out=list(c(cc0,paste(choose_gene,collapse = "+")))}else{ mat.out=list(rep(NA,2))}
      mat.out
    }
  
    
    library(doParallel) #并行处理包
    cl <- makeCluster(30)#makeCluster(detectCores())
    registerDoParallel(cl)
    mat=foreach(i=1:1000,.combine = "c")  %dopar% FGT3(i)
    mat.out3 <- data.frame(matrix(unlist(mat), nrow=length(mat), byrow=T),stringsAsFactors=FALSE)
    stopCluster(cl)  
  }
  mat.out3$number=sapply(stringr::str_split(mat.out3$X2,"\\+",simplify = F),function(x) length(x))
  mat.out3=unique(mat.out3)
  save(mat.out3,file="kaiqiaole_vfto.Rdata")#mat.outaa,
  load(file="kaiqiaole_vfto.Rdata")
  stringr::str_split(mat.out3[mat.out3$number==14,]$X2,"\\+",simplify=T)[1,]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # 
  # 
  # 
  # 
  # # rm(list = ls())  ## 魔幻操作，一键清空~ fpkm等于已经标准化完成了效果同TMM voom 这个是已经删除完了
  # library(WGCNA) #注意跑这个不要乱加载东西
  # library(stringr)
  # library(doParallel)
  # library(DT)
  # library(plyr)
  # library(gplots)
  # enableWGCNAThreads(nThreads = 8)
  # options(stringsAsFactors = F)
  # # load("C:/Users/zxh/Desktop/R/second/STAD_mrna_ENSG_fpkm.Rdata")
  # # load("C:/Users/zxh/Desktop/R/second/mRNA_fpkm_symbol_ebv.Rdata")
  # #load("rnatpmlog2001_Pancreas.Rdata") #原始350
  # #                                              load("rnatpmlog2001_Pancreas_pca.Rdata")
  # load("./rnatpm_Pancreas_WGCNA.Rdata")
  # 
  # 
  # rnaCounts=rnaCounts[,which(colnames(rnaCounts) %in% phe$sample)]#[mhighall_over,]#rownames(dePC)
  # rnatpm_PancreasT=t(rnaCounts)
  # rnatpm_PancreasT=as.data.frame(rnatpm_PancreasT)
  # 
  # dat=immune_timer2_CIBERSORT.ABS_dataset
  # # 每次都要检测数据
  # dat[1:4,1:4] 
  # dim(dat)
  # group_list=ifelse(rnatpm_PancreasT$ENSG00000140718>median(rnatpm_PancreasT$ENSG00000140718),'FTO_high','FTO_low')
  # table(group_list) #table函数，查看group_list中的分组个数
  # # 接下来走WGCNA流程： 
  # library(WGCNA)
  # ## WGCNA-step 1 : 构建输入数据，需要严格符合格式，否则会陷入调试代码的坑。
  # if(T){
  #   fpkm=dat
  #   fpkm[1:4,1:4]
  #   datTraits=data.frame(subtype=group_list)
  #   datTraits$gsm=colnames(dat)
  #   head(datTraits)
  #   table(datTraits$subtype)
  #   RNAseq_voom <- fpkm 
  #   ## 因为WGCNA针对的是基因进行聚类，而一般我们的聚类是针对样本用hclust即可，所以这个时候需要转置
  #   WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:22],])
  #   datExpr0 <- WGCNA_matrix  ## top 5000 mad genes
  #   datExpr <- datExpr0 
  #   ## 下面主要是为了防止临床表型与样本名字对不上
  #   sampleNames = rownames(datExpr);
  #   traitRows = match(sampleNames, datTraits$gsm)
  #   rownames(datTraits) = datTraits[traitRows, 'gsm']
  #   
  # }
  # #检验缺失的基因 无
  # if(T){
  #   gsg = goodSamplesGenes(datExpr, verbose = 3);
  #   gsg$allOK
  # }
  # if (!gsg$allOK){
  #   # Optionally, print the gene and sample names that were removed:
  #   if (sum(!gsg$goodGenes)>0)
  #     printFlush(paste("Removing genes:",
  #                      paste(colnames(datExpr)[!gsg$goodGenes], collapse = ",")));#names()
  #   if (sum(!gsg$goodSamples)>0)
  #     printFlush(paste("Removing samples:",
  #                      paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")));
  #   # Remove the offending genes and samples from the data:
  #   datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  # }
  # 
  # #检验离群样本 outliers 无
  # if(T){
  #   sampleTree = hclust(dist(datExpr), method = "average")
  #   pdf("kk.pdf",width = 80,height = 100)
  #   plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  #   dev.off()
  # }
  # #设置cutheight切之
  # if(F){
  #   clust = cutreeStatic(sampleTree, cutHeight = 400, minSize = 10)
  #   clust
  #   table(clust)
  #   #clust
  #   #0  1 
  #   #2 98 
  #   keepSamples = (clust==1)
  #   datExpr = datExpr[keepSamples, ]
  #   nGenes = ncol(datExpr)
  #   nSamples = nrow(datExpr)
  #   # save(datExpr, file = "FPKM-01-dataInput.RData")
  # }
  # #检验离群样本 outliers 牛逼的检测用z.k
  # if(T){
  #   suppressMessages(library(flashClust))
  #   # sample network based on squared Euclidean distance note that we
  #   # transpose the data
  #   A = adjacency(t(datExpr), type = "distance")
  #   # this calculates the whole network connectivity
  #   k = as.numeric(apply(A, 2, sum)) - 1
  #   # standardized connectivity
  #   Z.k = scale(k)
  #   # Designate samples as outlying if their Z.k value is below the threshold
  #   thresholdZ.k = -2.5  # often -2.5 可以用-0.4，这样删掉11/625个样本，但实际上没有必要，
  #   
  #   # the color vector indicates outlyingness (red)
  #   outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")
  #   table(outlierColor)
  #   # calculate the cluster tree using flahsClust or hclust
  #   suppressMessages(library(flashClust))
  #   sampleTree = flashClust(as.dist(1 - A), method = "average")
  #   # Convert traits to a color representation: where red indicates high
  #   # values
  #   datTraits_0=as.numeric(as.factor(datTraits$subtype))
  #   traitColors = data.frame(numbers2colors(datTraits_0, signed = FALSE))#原为datTraits
  #   dimnames(traitColors)[[2]] = paste(names(datTraits_0), "C", sep = "")
  #   datColors = data.frame(outlierC = outlierColor, traitColors)
  #   # Plot the sample dendrogram and the colors underneath.
  #   pdf("kk2.pdf",width = 80,height = 80)
  #   plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors, 
  #                       main = "Sample dendrogram and trait heatmap")
  #   dev.off()
  #   # plot of chunk unnamed-chunk-8
  #   # 
  #   # Caption: Cluster tree of mouse liver samples. The leaves of the tree correspond to the mice. The first color band underneath the tree indicates which arrays appear to be outlying. The second color band represents body weight (red indicates high value). Similarly, the remaining color-bands color-code the numeric values of physiologic traits.
  #   # 
  #   # The Figure shows the resulting cluster tree where each leaf corresponds to a mouse sample. The first color-band underneath the tree indicates which samples appear outlying (colored in red) according to a low value. The mouse labelled F2_221 is highly outlying (Z.k<-5), which is why we remove it from the subsequent analysis. The other color bands color-code physiological traits. Note that the outlying samples do not appear to have systematically different physiologic trait values. Although the outlying samples are suspicious, we will keep them in the subsequent analysis.
  #   # 
  #   # Remove outlying samples from expression and trait data
  #   remove.samples = Z.k < thresholdZ.k | is.na(Z.k)
  #   table(remove.samples)
  #   # the following 2 lines differ from what is written in the book
  #   datExpr = datExpr[!remove.samples, ]
  #   datTraits=as.data.frame(datTraits_0)
  #   datTraits = datTraits[!remove.samples, ]
  #   datTraits=data.frame(subtype=ifelse(datTraits==1,"FTO_low","FTO_high"))
  #   datTraits$gsm=rownames(datExpr)
  #   ## 下面主要是为了防止临床表型与样本名字对不上
  #   sampleNames = rownames(datExpr);
  #   traitRows = match(sampleNames, datTraits$gsm)
  #   rownames(datTraits) = datTraits[traitRows, 'gsm']
  # }
  # ## WGCNA-step 2 
  # ## 挑选最佳阈值，在  sft$powerEstimate
  # if(T){
  #   powers = c(1:10)#c(c(1:10), seq(from = 12, to=20, by=2))
  #   # Call the network topology analysis function
  #   sft = pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut=0.5, verbose = 5)
  #   #设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
  #   pdf("WGCNA-step2-beta-value.pdf",width = 8,height = 6)
  #   # Plot the results:
  #   ##sizeGrWindow(9, 5)
  #   par(mfrow = c(1,2));
  #   cex1 = 0.900;
  #   # Scale-free topology fit index as a function of the soft-thresholding power
  #   plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  #        xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
  #        main = paste("Scale independence"));
  #   text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  #        labels=powers,cex=cex1,col="red");
  #   # this line corresponds to using an R^2 cut-off of h
  #   abline(h=0.900,col="red")
  #   # Mean connectivity as a function of the soft-thresholding power
  #   plot(sft$fitIndices[,1], sft$fitIndices[,5],
  #        xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
  #        main = paste("Mean connectivity"))
  #   text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  #   dev.off()
  # }
  # sft$powerEstimate
  # sft[["fitIndices"]][["mean.k."]]
  # ## WGCNA-step3 构建加权共表达网络（Weight co-expression network)
  # ## 首先是一步法完成网络构建
  # if(T){
  #   net = blockwiseModules(
  #     datExpr,
  #     power = sft$powerEstimate,
  #     maxBlockSize = ncol(datExpr),
  #     TOMType = "unsigned", minModuleSize = 2,
  #     reassignThreshold = 0, mergeCutHeight = 0.5,
  #     numericLabels = TRUE, pamRespectsDendro = FALSE,
  #     saveTOMs = F,
  #     saveTOMFileBase = "AS-green-FPKM-TOM",
  #     verbose = 3
  #   )
  #   table(net$colors) 
  # }
  # 
  # ## 然后是分布法完成网络构建，仅供有探索精神的同学挑战。
  # 
  # ## 构建加权共表达网络分为两步：
  # ## 1. 计算邻近值，也是就是两个基因在不样品中表达量的表达相关系数(pearson correlation rho)，
  # ## 参考 2.b.2 in https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
  # ## 2. 计算topology overlap similarity (TOM)。 WGCNA认为，只通过计算两个基因的表达相关系数构建共表达网络是不足够的。
  # ## 于是他们用TOM表示两个基因在网络结构上的相似性，即两个基因如果具有相似的邻近基因，这两个基因更倾向于有相互作用。
  # ## 参考 2.b.3 in https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
  # 
  # 
  # if(F){
  #   #(1)网络构建 Co-expression similarity and adjacency 
  #   adjacency = adjacency(datExpr, power = sft$powerEstimate) 
  #   #(2) 邻近矩阵到拓扑矩阵的转换，Turn adjacency into topological overlap
  #   TOM = TOMsimilarity(adjacency);
  #   dissTOM = 1-TOM
  #   # (3) 聚类拓扑矩阵 Call the hierarchical clustering function
  #   geneTree = hclust(as.dist(dissTOM), method = "average");
  #   # Plot the resulting clustering tree (dendrogram)
  #   sizeGrWindow(12,9)
  #   ## 这个时候的geneTree与一步法的 net$dendrograms[[1]] 性质类似，但是还需要进行进一步处理
  #   plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
  #        labels = FALSE, hang = 0.04);
  #   #(4) 聚类分支的修整 dynamicTreeCut 
  #   # We like large modules, so we set the minimum module size relatively high:
  #   minModuleSize = 30;
  #   # Module identification using dynamic tree cut:
  #   dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
  #                               deepSplit = 2, pamRespectsDendro = FALSE,
  #                               minClusterSize = minModuleSize);
  #   table(dynamicMods)
  #   #4. 绘画结果展示
  #   # Convert numeric lables into colors
  #   dynamicColors = labels2colors(dynamicMods)
  #   table(dynamicColors)
  #   # Plot the dendrogram and colors underneath
  #   #sizeGrWindow(8,6)
  #   plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
  #                       dendroLabels = FALSE, hang = 0.03,
  #                       addGuide = TRUE, guideHang = 0.05,
  #                       main = "Gene dendrogram and module colors")
  #   #5. 聚类结果相似模块的融合，Merging of modules whose expression profiles are very similar
  #   #在聚类树中每一leaf是一个短线，代表一个基因，
  #   #不同分之间靠的越近表示有高的共表达基因，将共表达极其相似的modules进行融合
  #   # Calculate eigengenes
  #   MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  #   MEs = MEList$eigengenes
  #   # Calculate dissimilarity of module eigengenes
  #   MEDiss = 1-cor(MEs);
  #   # Cluster module eigengenes
  #   METree = hclust(as.dist(MEDiss), method = "average");
  #   # Plot the result
  #   #sizeGrWindow(7, 6)
  #   plot(METree, main = "Clustering of module eigengenes",
  #        xlab = "", sub = "")
  #   #选择有75%相关性的进行融合
  #   MEDissThres = 0.25
  #   # Plot the cut line into the dendrogram
  #   abline(h=MEDissThres, col = "red")
  #   # Call an automatic merging function
  #   merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  #   # The merged module colors
  #   mergedColors = merge$colors;
  #   # Eigengenes of the new merged modules:
  #   mergedMEs = merge$newMEs
  #   
  # }
  # 
  # ## WGCNA-step 4 ，看看模块颜色
  # if(T){
  #   
  #   # Convert labels to colors for plotting
  #   mergedColors = labels2colors(net$colors)
  #   table(mergedColors)
  #   moduleColors=mergedColors
  #   # Plot the dendrogram and the module colors underneath
  #   pdf("WGCNA-step4-genes-modules.pdf",width = 8,height = 6)
  #   plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
  #                       "Module colors",
  #                       dendroLabels = FALSE, hang = 0.03,
  #                       addGuide = TRUE, guideHang = 0.05)
  #   dev.off()
  #   ## assign all of the gene to their corresponding module 
  #   ## hclust for the genes.
  # }
  # # 如果感兴趣可以看看样本的形状分布情况。
  # # 一般来说，FTO_low和FTO_high肯定是泾渭分明的。
  # if(F){
  #   #明确样本数和基因
  #   nGenes = ncol(datExpr)
  #   nSamples = nrow(datExpr)
  #   #首先针对样本做个系统聚类
  #   datExpr_tree<-hclust(dist(datExpr), method = "average")
  #   par(mar = c(0,5,2,0))
  #   plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
  #        cex.axis = 1, cex.main = 1,cex.lab=1)
  #   ## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
  #   #针对前面构造的样品矩阵添加对应颜色
  #   sample_colors <- numbers2colors(as.numeric(factor(datTraits$subtype)), 
  #                                   colors = c( "red","green"),signed = FALSE)
  #   ## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目
  #   #  sample_colors <- numbers2colors( datTraits ,signed = FALSE)
  #   ## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，当然，这样给的颜色不明显，意义不大
  #   #10个样品的系统聚类树及性状热图
  #   par(mar = c(1,4,3,1),cex=0.8)
  #   
  #   pdf("sample-subtype-cluster.pdf",width = 8,height = 6)
  #   plotDendroAndColors(datExpr_tree, sample_colors,
  #                       groupLabels = colnames(sample),
  #                       cex.dendroLabels = 0.8,
  #                       marAll = c(1, 4, 3, 1),
  #                       cex.rowText = 0.01,
  #                       main = "Sample dendrogram and trait heatmap")
  #   dev.off()
  # }
  # # intramodularConnectivity
  # ## WGCNA-step 5 
  # ## 这一步主要是针对于连续变量，如果是分类变量，需要转换成连续变量方可使用
  # if(T){
  #   nGenes = ncol(datExpr)
  #   nSamples = nrow(datExpr)
  #   datTraits$subtype=as.factor(datTraits$subtype)
  #   design=model.matrix(~0+ datTraits$subtype)
  #   colnames(design)=levels(datTraits$subtype)
  #   design
  #   moduleColors <- labels2colors(net$colors)
  #   # Recalculate MEs with color labels
  #   MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  #   MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
  #   moduleTraitCor = cor(MEs, design , use = "p");
  #   moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  #   
  #   sizeGrWindow(10,6)
  #   # Will display correlations and their p-values
  #   textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
  #                      signif(moduleTraitPvalue, 1), ")", sep = "");
  #   dim(textMatrix) = dim(moduleTraitCor)
  #   pdf("WGCNA-step5-Module-trait-relationships.pdf",width = 8,height = 12)
  #   par(mar = c(6, 8.5, 3, 3));
  #   # Display the correlation values within a heatmap plot
  #   labeledHeatmap(Matrix = moduleTraitCor,
  #                  xLabels = colnames(design),
  #                  yLabels = names(MEs),
  #                  ySymbols = names(MEs),
  #                  colorLabels = FALSE,
  #                  colors = greenWhiteRed(50),
  #                  textMatrix = textMatrix,
  #                  setStdMargins = FALSE,
  #                  cex.text = 0.5,
  #                  zlim = c(-1,1),
  #                  main = paste("Module-trait relationships"))
  #   dev.off()
  # }
  # # 因为本例子是分类变量，而且就是泾渭分明的FTO_low和FTO_high，所以这个时候大家可以思考一下。
  # # 这个时候的WGCNA和差异分析的区别是什么。
  # 
  # ## WGCNA-step 6
  # 
  # if(T){
  #   # names (colors) of the modules
  #   modNames = substring(names(MEs), 3)
  #   geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  #   ## 算出每个模块跟基因的皮尔森相关系数矩
  #   ## MEs是每个模块在每个样本里面的
  #   ## datExpr是每个基因在每个样本的表达量
  #   MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  #   names(geneModuleMembership) = paste("MM", modNames, sep="");
  #   names(MMPvalue) = paste("p.MM", modNames, sep="");
  #   
  #   
  #   design
  #   ## 只有连续型性状才能只有计算
  #   ## 这里把是否是GC变量0,1进行数值化
  #   GC = as.data.frame(design[,2]);
  #   names(GC) = "GC"
  #   geneTraitSignificance = as.data.frame(cor(datExpr, GC, use = "p"));
  #   GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  #   names(geneTraitSignificance) = paste("GS.", names(GC), sep="");
  #   names(GSPvalue) = paste("p.GS.", names(GC), sep="");
  #   
  #   module = "turquoise"
  #   column = match(module, modNames);
  #   moduleGenes = moduleColors==module;
  #   pdf("WGCNA-step6-Module_membership-gene_significance.pdf",width = 8,height = 6)
  #   #sizeGrWindow(7, 7);
  #   par(mfrow = c(1,1));
  #   verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
  #                      abs(geneTraitSignificance[moduleGenes, 1]),
  #                      xlab = paste("Module Membership in", module, "module"),
  #                      ylab = "Gene significance for GC",
  #                      main = paste("Module membership vs. gene significance\n"),
  #                      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  #   abline(v=0.8,lwd=1,col="black")
  #   abline(h=0.2,lwd=1,col="black")
  #   dev.off()
  # }
  # 
  # ## WGCNA-step 7 ，对TOM矩阵进行可视化
  # if(T){
  #   nGenes = ncol(datExpr)
  #   nSamples = nrow(datExpr)
  #   geneTree = net$dendrograms[[1]]; 
  #   dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); 
  #   plotTOM = dissTOM^sft$powerEstimate; 
  #   diag(plotTOM) = NA; 
  #   #TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  #   nSelect = 1000#400
  #   # For reproducibility, we set the random seed
  #   set.seed(10);
  #   select = sample(nGenes, size = nSelect);
  #   selectTOM = dissTOM[select, select];
  #   # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  #   selectTree = hclust(as.dist(selectTOM), method = "average")
  #   selectColors = moduleColors[select];
  #   # Open a graphical window
  #   sizeGrWindow(9,9)
  #   # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
  #   # the color palette; setting the diagonal to NA also improves the clarity of the plot
  #   plotDiss = selectTOM^sft$powerEstimate;
  #   diag(plotDiss) = NA;
  #   
  #   pdf("WGCNA-step7-Network-heatmap.pdf",width = 8,height = 6)
  #   TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes",col=colorpanel(250,'red',"orange",'lemonchiffon'))
  #   dev.off()
  #   
  #   # Recalculate module eigengenes
  #   MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  #   ## 只有连续型性状才能只有计算
  #   ## 这里把是否属 GC 表型这个变量0,1进行数值化
  #   GC = as.data.frame(design[,2]);
  #   names(GC) = "GC"
  #   # Add the weight to existing module eigengenes
  #   MET = orderMEs(cbind(MEs, GC))
  #   # Plot the relationships among the eigengenes and the trait
  #   sizeGrWindow(5,7.5);
  #   
  #   par(cex = 0.9)
  #   pdf("WGCNA-step7-Eigengene-dendrogram.pdf",width = 8,height = 6)
  #   plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
  #                         = 90)
  #   dev.off()
  #   
  #   # Plot the dendrogram
  #   sizeGrWindow(6,6);
  #   par(cex = 1.0)
  #   ## 模块的进化树
  #   pdf("WGCNA-step7-Eigengene-dendrogram-hclust.pdf",width = 8,height = 6)
  #   plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
  #                         plotHeatmaps = FALSE)
  #   dev.off()
  #   # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  #   par(cex = 1.0)
  #   ## 性状与模块热
  #   
  #   pdf("WGCNA-step7-Eigengene-adjacency-heatmap.pdf",width = 8,height = 6)
  #   plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
  #                         plotDendrograms = FALSE, xLabelsAngle = 90)
  #   dev.off()
  #   
  # }
  # 
  # ## WGCNA-step 8；导出感兴趣模块的基因
  # if(T){
  #   # Select module
  #   module = "turquoise";
  #   # Select module probes
  #   probes = colnames(datExpr) ## 我们例子里面的probe就是基因
  #   head(moduleColors);head(probes)
  #   # 这个时候的 moduleColors 和  probes 是一一对应的。
  #   inModule = (moduleColors==module);
  #   modProbes = probes[inModule]; 
  #   head(modProbes)
  # }
  # 
  # ## WGCNA-step 9 ：导出指定的模块内部的基因之间的关联情况。
  # if(T){
  #   # Recalculate topological overlap
  #   TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); 
  #   # Select module
  #   module = "turquoise";
  #   # Select module probes
  #   probes = colnames(datExpr) ## 我们例子里面的probe就是基因
  #   inModule = (moduleColors==module);
  #   modProbes = probes[inModule]; 
  #   ## 也是提取指定模块的基因名
  #   # Select the corresponding Topological Overlap
  #   modTOM = TOM[inModule, inModule];
  #   dimnames(modTOM) = list(modProbes, modProbes)
  #   ## 模块对应的基因关系矩
  #   cyt = exportNetworkToCytoscape(
  #     modTOM,
  #     edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  #     nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  #     weighted = TRUE,
  #     threshold = 0.02,
  #     nodeNames = modProbes, 
  #     nodeAttr = moduleColors[inModule]
  #   );
  #   # 导出指定的模块内部的基因之间的关联情况。
  #   head(cyt$edgeData)
  #   head(cyt$nodeData)
  # }
  # ## 对模块的基因集批量GO/KEGG富集分析
  # 
  # FilterGenes_spe <- ((abs(geneModuleMembership[moduleGenes, column]) > 0.8)
  #                     & (abs(geneTraitSignificance[moduleGenes, 1]) > 0.2) )
  # table(FilterGenes_spe)
  # trait_hubGenes_spe <- modProbes[FilterGenes_spe]
  # hubGene=as.data.frame(trait_hubGenes_spe)
  # if(length(hubGene)!=0){write.csv(trait_hubGenes_spe,file='trait_hubGenes_spe.csv',quote = T)
  # }else{
  #   print("no hub under 0.8 MM 0.2 GS")}
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # # tcga_only_expr_mean=t(scale(t(exprSet_quant_with_clinical)))
  # # 
  # # ui=cbind(tcga_only_expr_mean[mhighall_over,],
  # #          genes_expr_mean_EMTAB6134[mhighall_over,],
  # #          ICGC_PACA_CA_seq_log_filter_dornor_mean[mhighall_over,],
  # #          ICGC_PACA_AU_array_norm_dornor_mean[mhighall_over,],
  # #          genes_expr_mean_GSE71729_rna[mhighall_over,],
  # #          genes_expr_mean_GSE21501_rna[mhighall_over,]
  # # )