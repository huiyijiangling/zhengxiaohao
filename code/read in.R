# read in 
# expr <- read.table("LIHC.TPM.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)






#####

library(tidyverse) # 用于读取MAF文件
library(ISOpureR) # 用于纯化表达谱
library(impute) # 用于KNN填补药敏数据
library(pRRophetic) # 用于药敏预测
library(SimDesign) # 用于禁止药敏预测过程输出的信息
library(ggplot2) # 绘图
library(cowplot) # 合并图像
library(oncoPredict) # 合并图像
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}  


# # 输入文件
# 
# 文件较大，已上传至微云，请点击链接下载<https://share.weiyun.com/c9oY6n0T>
#   
# LIHC.TPM.txt，基因表达矩阵。计算PPS得分和药敏预测，都基于这个表达矩阵。实际应用时把这个表达矩阵替换成你自己的数据。
# 
# data_mutations_mskcc.txt，maf突变数据，用于提取癌症亚型，例文提取的是TP53-mutant patients。


# 1.读取肝癌TPM表达谱
load(r"(C:\Users\zxh\Desktop\R\meta collect\ensg_symbol_relation.Rdata)")
rownames(ensg_symbol_trans)=stringr::str_split(rownames(ensg_symbol_trans),"\\.",simplify=T)[,1]
load(r"(C:\Users\zxh\Desktop\R\meta collect\tcga_pcawg_rnaseq.Rdata)")
expr=TCGA_tpm_symbol[["PAAD"]]

before_multiple_surviaval <- function(genes, rna.expr, metaMatrix){
  library(survminer)
  library(survival)
  samples = intersect(colnames(rna.expr), metaMatrix$sample)
  exprDa=rna.expr[genes,samples]
  clinicalDa=metaMatrix[match(samples,metaMatrix$sample),]
  colnames(exprDa)==rownames(clinicalDa)
  dat=cbind(t(exprDa),clinicalDa)
  return(dat)
}
dat=before_multiple_surviaval("ENSG00000140718",TCGA_tpm_ensg[["PAAD"]], metaMatrix=pcawg_clinical_CDR)
# duplicated(dat$sample)
# dat=merge(dat,PAAD_knownclass_TCGA,by="sample")
dat$FTO=ifelse(dat$ENSG00000140718>median(dat$ENSG00000140718,na.rm = T),"FTO-high","FTO-low")
##########################
load(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\ICGC_PACA_CA_metaMatrix_duct_surv.Rdata)")
expr=ICGC_PACA_CA_seq_log_filter_dornor
# 加载基因注释文件，用于基因ID转换
comgene <- intersect(rownames(expr),rownames(ensg_symbol_trans))
expr <- expr[comgene,]
expr$gene <- ensg_symbol_trans[comgene,"gene_name"];
expr <- expr[!duplicated(expr$gene),];
rownames(expr) <- expr$gene; 
expr <- expr[,-ncol(expr)]
##########################
load(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\ICGC_PACA_AU_metaMatrix_duct_surv.Rdata)")
expr=ICGC_PACA_AU_array_norm_dornor
# 加载基因注释文件，用于基因ID转换
comgene <- intersect(rownames(expr),rownames(ensg_symbol_trans))
expr <- expr[comgene,]
expr$gene <- ensg_symbol_trans[comgene,"gene_name"];
expr <- expr[!duplicated(expr$gene),];
rownames(expr) <- expr$gene; 
expr <- expr[,-ncol(expr)]
###########
if(F){

normsam <- colnames(expr[,which(substr(colnames(expr),14,15) == "11")])
tumosam <- colnames(expr[,which(substr(colnames(expr),14,15) == "01")])

# 2.读取maf突变文件(于cBioPortal下载)
# maf <- read_tsv("data_mutations_mskcc.txt", comment = "#")
# maf$Tumor_Sample_Barcode <- paste0("LIHC",substr(maf$Tumor_Sample_Barcode,8,15))
if(F){
maf <- read_tsv(r"(C:\Users\zxh\Desktop\R\meta collect\pcawg_publish\source\mc3.v0.2.8.PUBLIC.maf.gz)", comment = "#")

# 提取既有表达数据又有突变数据的肿瘤样本
maf$Tumor_Sample_Barcode=substr(maf$Tumor_Sample_Barcode,1,15)
tumosam <- intersect(tumosam,unique(maf$Tumor_Sample_Barcode)) 
maf <- maf[which(maf$Tumor_Sample_Barcode %in% tumosam),]
expr <- expr[,c(tumosam,normsam)]

# 3.提取TP53突变信息，并创建样本注释
tp53 <- c()
for (i in tumosam) {
  tmp <- maf[which(maf$Tumor_Sample_Barcode == i),]
  if(is.element("TP53", tmp$Hugo_Symbol)) { # 如果存在TP53
    tp53 <- c(tp53,1) # 记录1
  } else {
    tp53 <- c(tp53,0) # 否则记录0
  }
}
names(tp53) <- tumosam

# 取出有TP53突变的患者
tp53.mutsam <- names(tp53[which(tp53 == 1)]) 

signature <- read.table("17gene.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

pps <- as.numeric(apply(t(log2(expr[rownames(signature),tp53.mutsam] + 1)), 1, function(x) {x %*% signature$Coefficient}))

# 2.标准化，把pps处理到0-1之间
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
npps <- range01(pps)

# 3.创建样本信息
Sinfo <- data.frame(PPS = npps,
                    TP53 = tp53[tp53.mutsam],
                    row.names = tp53.mutsam,
                    stringsAsFactors = F)
head(Sinfo)

# 把pps保存到文件
write.csv(Sinfo, "output_PPS.csv", quote = F)
}

# 如果想运行就把`runpure <- F`改为`runpure <- T`

normexpr <- as.matrix(expr[,normsam])
tumoexpr <- as.matrix(expr[,tumosam])#tp53.mutsam
}
tumoexpr=expr
runpure <- F # 如果想运行就把这个改为T
if(runpure) {
  set.seed(123)
  # Run ISOpureR Step 1 - Cancer Profile Estimation
  ISOpureS1model <- ISOpure.step1.CPE(tumoexpr, normexpr)
  # For reproducible results, set the random seed
  set.seed(456);
  # Run ISOpureR Step 2 - Patient Profile Estimation
  ISOpureS2model <- ISOpure.step2.PPE(tumoexpr,normexpr,ISOpureS1model)
  pure.tumoexpr <- ISOpureS2model$cc_cancerprofiles
}

if(!runpure) {
  pure.tumoexpr <- tumoexpr
}

# 1.制作CTRP AUC矩阵，保存到CTRP_AUC.txt文件
auc <- read.table("CTRP_AUC_raw.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 3
auc$comb <- paste(auc$master_cpd_id,auc$master_ccl_id,sep = "-")
auc <- apply(auc[,"area_under_curve",drop = F], 2, function(x) tapply(x, INDEX=factor(auc$comb), FUN=max, na.rm=TRUE)) # 重复项取最大AUC
auc <- as.data.frame(auc)
auc$master_cpd_id <- sapply(strsplit(rownames(auc),"-",fixed = T),"[",1)
auc$master_ccl_id <- sapply(strsplit(rownames(auc),"-",fixed = T),"[",2)
auc <- reshape(auc, 
               direction = "wide",
               timevar = "master_cpd_id",
               idvar = "master_ccl_id")
colnames(auc) <- gsub("area_under_curve.","",colnames(auc),fixed = T)
ctrp.ccl.anno <- read.table("CTRP_ccl_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 1
ctrp.cpd.anno <- read.delim("CTRP_cpd_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 2
# 保存到文件
write.table(auc,"CTRP_AUC.txt",sep = "\t",row.names = F,col.names = T,quote = F)

# 2.加载药敏AUC矩阵并进行数据处理
ctrp.auc <- read.table("CTRP_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
prism.auc <- read.delim("PRISM_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # 数据来自https://depmap.org/portal/download/ Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4
prism.ccl.anno <- prism.auc[,1:5] # 前5列为细胞系注释信息
prism.auc <- prism.auc[,-c(1:5)]
# ctrp.auc 横的是细胞 列的是药
# prism.auc 横的是细胞 列的是药


#理论上这个可加上了
## a. 移除缺失值大于20%的药物
ctrp.auc <- ctrp.auc[,apply(ctrp.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(ctrp.auc)]
ctrp.auc <- ctrp.auc[apply(ctrp.auc,1,function(x) sum(is.na(x))) < 0.5*ncol(ctrp.auc),]

prism.auc <- prism.auc[,apply(prism.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(prism.auc)]
prism.auc <- prism.auc[apply(prism.auc,1,function(x) sum(is.na(x))) < 0.5*ncol(prism.auc),]

## b. 移除CTRP数据里源自haematopoietic_and_lymphoid_tissue的细胞系
if(F){
rmccl <- paste0("CCL",na.omit(ctrp.ccl.anno[which(ctrp.ccl.anno$ccle_primary_site == "haematopoietic_and_lymphoid_tissue"),"master_ccl_id"]))
}
rownames(ctrp.auc) <- paste0("CCL",rownames(ctrp.auc))
if(F){
ctrp.auc <- ctrp.auc[setdiff(rownames(ctrp.auc),rmccl),]
}
if(T){
## c. KNN填补缺失值
ctrp.auc.knn <- impute.knn(as.matrix(ctrp.auc))$data
prism.auc.knn <- impute.knn(as.matrix(prism.auc))$data
}
# ctrp.auc.knn=ctrp.auc
# prism.auc.knn=prism.auc
## d. 数据量级修正（与作者沟通得知） 0是最敏感1是最不敏感1，1-auc是敏感度
ctrp.auc.knn <- ctrp.auc.knn/ceiling(max(ctrp.auc.knn)) # 参考Expression Levels of Therapeutic Targets as Indicators of Sensitivity to Targeted Therapeutics (2019, Molecular Cancer Therapeutics)
prism.auc.knn <- prism.auc.knn/ceiling(max(prism.auc.knn))

# 药敏预测

# 加载CCLE细胞系的表达谱，作为训练集
ccl.expr <- read.table("CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 


# 把基因的ensembl ID转换为gene symbol
ccl.expr <- ccl.expr[,-1]
rownames(ccl.expr) <- sapply(strsplit(rownames(ccl.expr),".",fixed = T),"[",1)

# 加载基因注释文件，用于基因ID转换
Ginfo <- read.table("overlapTable27.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # 参考FigureYa34count2FPKMv2制作的基因注释文件
comgene <- intersect(rownames(ccl.expr),rownames(Ginfo))
ccl.expr <- ccl.expr[comgene,]
ccl.expr$gene <- Ginfo[comgene,"genename"];
ccl.expr <- ccl.expr[!duplicated(ccl.expr$gene),];
rownames(ccl.expr) <- ccl.expr$gene; 
ccl.expr <- ccl.expr[,-ncol(ccl.expr)]


# 下面用pRRophetic包里的calcPhenotype函数，分别基于CTRP和PRISM，估计每个样本的drug response。

## CTRP

keepgene <- apply(ccl.expr, 1, mad) > 0.5 # 保留表达值有效的基因
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1) # 重置细胞系名
trainPtype <- as.data.frame(ctrp.auc.knn)
ccl.name <- ccl.miss <- c() # 替换细胞系名
for (i in rownames(trainPtype)) {
  if(!is.element(gsub("CCL","",i),ctrp.ccl.anno$master_ccl_id)) {
    cat(i,"\n")
    ccl.miss <- c(ccl.miss, i) # 没有匹配到的细胞系
    ccl.name <- c(ccl.name, i) # 插入未匹配的细胞系
  } else {
    ccl.name <- c(ccl.name,  ctrp.ccl.anno[which(ctrp.ccl.anno$master_ccl_id == gsub("CCL","",i)),"ccl_name"]) # 插入匹配的细胞系
  }
}

cpd.name <- cpd.miss <- c() # 替换药物名
for (i in colnames(trainPtype)) {
  if(!is.element(i,ctrp.cpd.anno$master_cpd_id)) {
    cat(i,"\n")
    cpd.miss <- c(cpd.miss, i) # 没有匹配到的药物
    cpd.name <- c(cpd.name, i) # 插入未匹配的药物
  } else {
    cpd.name <- c(cpd.name,  ctrp.cpd.anno[which(ctrp.cpd.anno$master_cpd_id == i),"cpd_name"]) # 插入匹配的药物
  }
}

rownames(trainPtype) <- ccl.name
trainPtype <- trainPtype[setdiff(rownames(trainPtype),ccl.miss),] # 去除未匹配的细胞系
colnames(trainPtype) <- cpd.name
trainPtype <- trainPtype[,setdiff(colnames(trainPtype),cpd.miss)] # 去除未匹配的药物
comccl <- intersect(rownames(trainPtype),colnames(trainExpr)) # 提取有表达且有药敏的细胞系
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

# 测试集
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5 # 纯化的测试集取表达稳定的基因
testExpr <- log2(pure.tumoexpr[keepgene,] + 1) # 表达谱对数化
# 取训练集和测试集共有的基因
comgene <- intersect(rownames(trainExpr),rownames(testExpr)) 
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]


# 岭回归预测药物敏感性

outTab <- NULL
# 循环很慢，请耐心
#原文算到这里我开始改了
if(F){
  for (i in 1:ncol(trainPtype)) { 
    display.progress(index = i,totalN = ncol(trainPtype))
    d <- colnames(trainPtype)[i]
    tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于CTRP的AUC可能有0值，因此加一个较小的数值防止报错
    
    # 岭回归预测药物敏感性
    ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                    trainingPtype = tmp,
                                    testExprData = testExpr,
                                    powerTransformPhenotype = F,
                                    selection = 1))
    ptypeOut <- 2^ptypeOut - 0.00001 # 反对数
    outTab <- rbind.data.frame(outTab,ptypeOut)
  }
}
ptypeOut <- oncoPredict::calcPhenotype(trainingExprData = trainExpr,
                                       trainingPtype = as.matrix(trainPtype),
                                       testExprData = as.matrix(testExpr),
                                       removeLowVaryingGenes = 0.2, 
                                       minNumSamples = 10,
                                       batchCorrect = "eb", ##batchCorrect options: "eb" for ComBat, "qn" for quantiles normalization, "standardize", or "none"
                                       #"eb" is good to use when you use microarray training data to build models on microarray testing data.
                                       #"standardize is good to use when you use microarray training data to build models on RNA-seq testing data (this is what Paul used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data, see methods section of that paper for rationale)
                                       removeLowVaringGenesFrom = "homogenizeData",##Options are 'homogenizeData' and 'rawData'
                                       #homogenizeData is likely better if there is ComBat batch correction,
                                       printOutput = TRUE, 
                                       powerTransformPhenotype = F,
                                       selection = 1,
                                       #Indicate whether or not you'd like to use PCA for feature/gene reduction. Options are 'TRUE' and 'FALSE'.
                                       #Note: If you indicate 'report_pca=TRUE' you need to also indicate 'pca=TRUE'
                                       pcr = FALSE,
                                       
                                       #Indicate whether you want to output the principal components. Options are 'TRUE' and 'FALSE'.
                                       report_pc = FALSE,
                                       
                                       #Indicate if you want correlation coefficients for biomarker discovery. These are the correlations between a given gene of interest across all samples vs. a given drug response across samples.
                                       #These correlations can be ranked to obtain a ranked correlation to determine highly correlated drug-gene associations.
                                       cc = T,
                                       
                                       #Indicate whether or not you want to output the R^2 values for the data you train on from true and predicted values.
                                       #These values represent the percentage in which the optimal model accounts for the variance in the training data.
                                       #Options are 'TRUE' and 'FALSE'.
                                       rsq = T,
                                       
                                       #Indicate percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
                                       percent = 80)
# outTab=ptypeOut


outTab=read.csv(r"(C:\Users\zxh\Desktop\R\paad-tcga-gtex\samlltu\FigureYa212drugTargetV2\sum\TGCA\CTRP\calcPhenotype_Output\DrugPredictions.csv)",row.names = 1)
ctrp.pred.auc <- outTab
# # original<-getwd()
# # wd<-tempdir()
# # savedir<-setwd(wd)
# 
outTab=t(outTab)
outTab=as.data.frame(outTab)
# dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc <- outTab
rownames(ctrp.pred.auc)=tolower(rownames(ctrp.pred.auc))
## PRISM

# 数据准备


keepgene <- apply(ccl.expr, 1, mad) > 0.5
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1)
trainPtype <- as.data.frame(prism.auc.knn)
rownames(trainPtype) <- prism.ccl.anno[rownames(trainPtype),"cell_line_display_name"]
#colnames(trainPtype) <- sapply(strsplit(colnames(trainPtype)," (",fixed = T), "[",1)
comccl <- intersect(rownames(trainPtype),colnames(trainExpr))
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

# 测试集
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5
testExpr <- log2(pure.tumoexpr[keepgene,] + 1)
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]



outTab <- NULL
# 循环很慢，请耐心
# for (i in 1:ncol(trainPtype)) { 
#   display.progress(index = i,totalN = ncol(trainPtype))
#   d <- colnames(trainPtype)[i]
#   tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于PRISM的AUC可能有0值，因此加一个较小的数值防止报错
#   ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
#                                   trainingPtype = tmp,
#                                   testExprData = testExpr,
#                                   powerTransformPhenotype = F,#Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#                                   selection = 1))
#   ptypeOut <- 2^ptypeOut - 0.00001 # 反对数
#   outTab <- rbind.data.frame(outTab,ptypeOut)
# }
ptypeOut <- oncoPredict::calcPhenotype(trainingExprData = trainExpr,
                                       trainingPtype = as.matrix(trainPtype),
                                       testExprData = as.matrix(testExpr),
                                       removeLowVaryingGenes = 0.2, 
                                       minNumSamples = 10,
                                       batchCorrect = "eb", ##batchCorrect options: "eb" for ComBat, "qn" for quantiles normalization, "standardize", or "none"
                                       #"eb" is good to use when you use microarray training data to build models on microarray testing data.
                                       #"standardize is good to use when you use microarray training data to build models on RNA-seq testing data (this is what Paul used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data, see methods section of that paper for rationale)
                                       removeLowVaringGenesFrom = "homogenizeData",##Options are 'homogenizeData' and 'rawData'
                                       #homogenizeData is likely better if there is ComBat batch correction,
                                       printOutput = TRUE, 
                                       powerTransformPhenotype = F,
                                       selection = 1,
                                       #Indicate whether or not you'd like to use PCA for feature/gene reduction. Options are 'TRUE' and 'FALSE'.
                                       #Note: If you indicate 'report_pca=TRUE' you need to also indicate 'pca=TRUE'
                                       pcr = FALSE,
                                       
                                       #Indicate whether you want to output the principal components. Options are 'TRUE' and 'FALSE'.
                                       report_pc = FALSE,
                                       
                                       #Indicate if you want correlation coefficients for biomarker discovery. These are the correlations between a given gene of interest across all samples vs. a given drug response across samples.
                                       #These correlations can be ranked to obtain a ranked correlation to determine highly correlated drug-gene associations.
                                       cc = T,
                                       
                                       #Indicate whether or not you want to output the R^2 values for the data you train on from true and predicted values.
                                       #These values represent the percentage in which the optimal model accounts for the variance in the training data.
                                       #Options are 'TRUE' and 'FALSE'.
                                       rsq = T,
                                       
                                       #Indicate percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
                                       percent = 80)


outTab=read.csv(r"(C:\Users\zxh\Desktop\R\paad-tcga-gtex\samlltu\FigureYa212drugTargetV2\sum\TGCA\PRISM\calcPhenotype_Output\DrugPredictions.csv)",row.names = 1)
prism.pred.auc <- outTab
outTab=t(outTab)
outTab=as.data.frame(outTab)
prism.pred.auc <- outTab
prism.pred.auc=prism.pred.auc[!duplicated(stringr::str_split(rownames(prism.pred.auc),"\\.\\.",simplify = T)[,1]),]
rownames(prism.pred.auc)=stringr::str_split(rownames(prism.pred.auc),"\\.\\.",simplify = T)[,1]
rownames(prism.pred.auc)=tolower(rownames(prism.pred.auc))
# 确定潜在药物靶标

Sinfo=as.data.frame(t(testExpr))
Sinfo=Sinfo[intersect(names(ctrp.pred.auc),names(prism.pred.auc)),]
# Sinfo=as.data.frame(testExpr["FTO",])
rownames(ctrp.pred.auc)%in%rownames(prism.pred.auc)

save(ctrp.pred.auc,prism.pred.auc,Sinfo,file="TCGA_pred.Rdata")
load(file="TCGA_pred.Rdata")
#只要一样的
common_drug=intersect(rownames(ctrp.pred.auc),rownames(prism.pred.auc))
ctrp.pred.auc=ctrp.pred.auc[common_drug,]
prism.pred.auc=prism.pred.auc[common_drug,]

write.csv(ctrp.pred.auc,file ="ctrp.pred.auc.csv" )
write.csv(prism.pred.auc,file ="prism.pred.auc.csv" )
Sinfo$PPS=Sinfo$FTO
top.pps <- Sinfo[Sinfo$PPS >= quantile(Sinfo$PPS,probs = seq(0,1,0.25))[4],] # 定义上十分位的样本
bot.pps <- Sinfo[Sinfo$PPS <= quantile(Sinfo$PPS,probs = seq(0,1,0.25))[2],] # 定义下十分位的样本
rownames(Sinfo)%in%colnames(ctrp.pred.auc)


## 1.差异药敏分析
ctrp.log2fc <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  display.progress(index = i,totalN = nrow(ctrp.pred.auc))
  d <- rownames(ctrp.pred.auc)[i]
  a <- mean(as.numeric(ctrp.pred.auc[d,rownames(top.pps)])) # 上十分位数的AUC均值
  b <- mean(as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])) # 下十分位数的AUC均值
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  ctrp.log2fc <- c(ctrp.log2fc,log2fc)
}
candidate.ctrp <- ctrp.log2fc#[ctrp.log2fc > 0.2] # 这里我调整了阈值，控制结果数目

prism.log2fc <- c()
for (i in 1:nrow(prism.pred.auc)) {
  display.progress(index = i,totalN = nrow(prism.pred.auc))
  d <- rownames(prism.pred.auc)[i]
  a <- mean(as.numeric(prism.pred.auc[d,rownames(top.pps)])) # 上十分位数的AUC均值
  b <- mean(as.numeric(prism.pred.auc[d,rownames(bot.pps)])) # 下十分位数的AUC均值
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  prism.log2fc <- c(prism.log2fc,log2fc)
}
candidate.prism <- prism.log2fc#[prism.log2fc > 0.2] # 这里我调整了阈值，控制结果数目

## 2.Spearman相关性分析，用于绘制左图

ctrp.cor <- ctrp.cor.p <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  display.progress(index = i,totalN = nrow(ctrp.pred.auc))
  d <- rownames(ctrp.pred.auc)[i]
  a <- as.numeric(ctrp.pred.auc[d,rownames(Sinfo)]) 
  b <- as.numeric(Sinfo$PPS)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  ctrp.cor <- c(ctrp.cor,r)
  ctrp.cor.p <- c(ctrp.cor.p,p)
}
candidate.ctrp2 <- ctrp.cor#[ctrp.cor < -0.4]  # 这里我调整了阈值，控制结果数目
ctrp.candidate <- intersect(names(candidate.ctrp),names(candidate.ctrp2))

prism.cor <- prism.cor.p <- c()
for (i in 1:nrow(prism.pred.auc)) {
  display.progress(index = i,totalN = nrow(prism.pred.auc))
  d <- rownames(prism.pred.auc)[i]
  a <- as.numeric(prism.pred.auc[d,rownames(Sinfo)]) 
  b <- as.numeric(Sinfo$PPS)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  prism.cor <- c(prism.cor,r)
  prism.cor.p <- c(prism.cor.p,p)
}
candidate.prism2 <- prism.cor#[prism.cor < -0.35]  
prism.candidate <- intersect(names(candidate.prism),names(candidate.prism2))
```

# 开始画图

## 1. 左侧相关性图

```{r}
# 设置颜色
darkblue <- "#0772B9"
lightblue <- "#48C8EF"

cor.data <- data.frame(drug = ctrp.candidate,
                       r = ctrp.cor[ctrp.candidate],
                       p = -log10(ctrp.cor.p[ctrp.candidate]+10^(-10)))
p1 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.3, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank())

cor.data <- data.frame(drug = prism.candidate,
                       r = prism.cor[prism.candidate],
                       p = -log10(prism.cor.p[prism.candidate]))
cor.data$drug <- sapply(strsplit(cor.data$drug," (",fixed = T), "[",1)
dplyr::top_n(cor.data,10,wt=cor.data$r)
p2 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.3, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank())


## 2.右侧箱型图


ctrp.boxdata <- NULL
for (d in ctrp.candidate) {
  a <- as.numeric(ctrp.pred.auc[d,rownames(top.pps)]) 
  b <- as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  ctrp.boxdata <- rbind.data.frame(ctrp.boxdata,
                                   data.frame(drug = d,
                                              auc = c(a,b),
                                              p = p,
                                              s = s,
                                              group = rep(c("High PPS","Low PPS"),c(nrow(top.pps),nrow(bot.pps))),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
p3 <- ggplot(ctrp.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=ctrp.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
            label=ctrp.boxdata$s,
            data=ctrp.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(darkblue, lightblue)) + 
  scale_color_manual(values = c(darkblue, lightblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) 
dat <- ggplot_build(p3)$data[[1]]

p3 <- p3 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)

prism.boxdata <- NULL
for (d in prism.candidate) {
  a <- as.numeric(prism.pred.auc[d,rownames(top.pps)]) 
  b <- as.numeric(prism.pred.auc[d,rownames(bot.pps)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  prism.boxdata <- rbind.data.frame(prism.boxdata,
                                    data.frame(drug = d,
                                               auc = c(a,b),
                                               p = p,
                                               s = s,
                                               group = rep(c("High PPS","Low PPS"),c(nrow(top.pps),nrow(bot.pps))),
                                               stringsAsFactors = F),
                                    stringsAsFactors = F)
}
prism.boxdata$drug <- sapply(strsplit(prism.boxdata$drug," (",fixed = T), "[",1)

p4 <- ggplot(prism.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=prism.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
            label=prism.boxdata$s,
            data=prism.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(darkblue, lightblue)) + 
  scale_color_manual(values = c(darkblue, lightblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 10),
        legend.position = "bottom",
        legend.title = element_blank())
dat <- ggplot_build(p4)$data[[1]]

p4 <- p4 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)


## 3. 合并图像
# p2, p4,
```{r, fig.width=8, fig.height=8}
ppp=plot_grid(p1,p2,p3,p4,  labels=c("A", "", "B", ""), 
          ncol=2, 
          rel_widths = c(2, 2)) #左右两列的宽度比例
save(ppp,file = "ppp.Rdata")
ggsave(filename = "drug target1.pdf",width = 40,height = 40)

# 保存镜像
#save.image("drug.RData")
```

# Session Info

```{r}
sessionInfo()
```