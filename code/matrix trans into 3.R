rm(list=ls())
gc()
options(stringsAsFactors = F)
library(Seurat)
library("DropletUtils")
# GSM4730268 ok 要反一下
setwd(r"(C:\Users\zxh\Desktop\R\scCancer\GSE156405\GSM4730265\filtered_feature_bc_matrix)")
R.utils::gzip(file = "barcodes.tsv",remove=TRUE,overwrite=T)
fs=list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE156405\GSM4730268)",pattern = '.csv',full.names = T)
fs
ensg=list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE156405\)",pattern = 'features.tsv.gz',full.names = T,recursive = T)
ensg
kk=list()
kk=lapply(ensg,function(x)readr::read_tsv(x,col_names = F))
kk=data.table::rbindlist(kk)
kk=unique(kk)
# GSM4730268 ok
fs=list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE156405\GSM4730268)",pattern = '.csv',full.names = T)
fs
# GSE165399 ok xinalinai
fs=list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE165399)",pattern = 'matrix.txt.gz',full.names = T)
# GSE137766 ok
fs=list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE137766)",pattern = 'umi.csv.gz',full.names = T)
fs
# GSE134355 ok 只处理了 GSM4008637
fs=list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE134355)",pattern = 'GSM4008637_Adult-Pancreas1_dge.txt',full.names = T)
fs
# GSE84133 ok 只处理了 GSM2230760 GSM2230759 GSM2230758 GSM2230757
fs=list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE84133)",pattern = '_umifm_counts.csv.gz',full.names = T)
fs=fs[1:4]
fs
#
# GSE156405 ok
fs=list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE156405)",pattern = '_filtered_feature_bc_matrix',full.names = T)
fs
list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE156405)",pattern = 'GSM')
# GSE156405 ok
fs=list.files(r"(H:\download.big.ac.cn\gsa\CRA001160)",pattern = 'count-matrix.txt',full.names = T)
fs
# tcga
load(r"(C:\Users\zxh\Desktop\R\meta collect\ensg_symbol_relation.Rdata)")
rownames(ensg_symbol_trans)=stringr::str_split(rownames(ensg_symbol_trans),"\\.",simplify=T)[,1]
load(r"(C:\Users\zxh\Desktop\R\meta collect\tcga_pcawg_rnaseq_counts.Rdata)")

fs=TCGA_counts_symbol[["PAAD"]]

# rownames(TCGA_counts_symbol[[1]])==

ensg_symbol_trans=ensg_symbol_trans[!duplicated(ensg_symbol_trans$gene_name),]
trans=data.frame(rownames(TCGA_counts_symbol[[1]]),rownames(TCGA_counts_symbol[[1]]))
trans=merge(ensg_symbol_trans,trans,by.x="gene_name",by.y="rownames.TCGA_counts_symbol..1....1")
trans[,c(1,3)]=NULL
trans$gene_type="Gene Expression"
rownames(TCGA_counts_symbol[[1]])==trans$rownames.TCGA_counts_symbol..1...
trans=trans[order(match(trans$rownames.TCGA_counts_symbol..1...,rownames(TCGA_counts_symbol[[1]]))),]
rownames(TCGA_counts_symbol[[1]])==trans$rownames.TCGA_counts_symbol..1...
x=1
#处理TCGA
sceList_tcga <-  lapply(1:33, function(x){
  p=names(TCGA_counts_symbol)[x]
  print(p)
  ct=as.matrix(TCGA_counts_symbol)[[x]]
  dir.create(file.path(p))
  dir.create(file.path(paste0(p,'/filtered_feature_bc_matrix/')))#filtered_feature_bc_matrix 
  data.table::fwrite(trans,file = paste0(p,'/filtered_feature_bc_matrix/features.tsv.gz'),
                     quote = F,sep = '\t',
                     col.names = F,row.names = F)
  #genes.tsv = symbol,symbol,GE, features= ensg,symbol,GE
  data.table::fwrite(data.frame(colnames(ct)),file = paste0(p,'/filtered_feature_bc_matrix/barcodes.tsv.gz'),quote = F,sep = '\t',
                     col.names = F,row.names = F)
  file=paste0(p,'/filtered_feature_bc_matrix/matrix.mtx')
  sink(file)
  cat("%%MatrixMarket matrix coordinate integer general\n")
  cat("%\n")
  cat(paste(nrow(ct),ncol(ct),sum(ct>0),"\n")) 
  sink()
  
  #再写入表达量信息
  # tmp1=do.call(rbind,lapply(1:ncol(ct),function(i){
  #   return(data.frame(row=1:nrow(ct),
  #                     col=i,
  #                     exp=ct[,i]))
  # }) )
  tmp=data.table::rbindlist(lapply(1:ncol(ct),function(i){
    return(data.frame(row=1:nrow(ct),
                      col=i,
                      exp=ct[,i]))
  }))
  tmp=tmp[tmp$exp>0,]
  head(tmp)
  data.table::fwrite(tmp,file = paste0(p,'/filtered_feature_bc_matrix/matrix.mtx'),quote = F,sep = '\t',
                     col.names = F,row.names = F,append = T)
  R.utils::gzip(file = paste0(p,'/filtered_feature_bc_matrix/matrix.mtx'),remove=TRUE,overwrite=T)
})

#处理文本型
sceList <-  lapply(fs, function(x){
  x=fs
  raw.data=data.table::fread(x)
  raw.data=as.data.frame(raw.data)
  raw.data[1:4,1:5]
  #一般
  if(T){
  rownames(raw.data)=raw.data[,1]
  raw.data=raw.data[,-1]
  }
  #万一barcode和gene是反的
  if(F){
  raw.data=t(raw.data)  
  }
  if(F){
  metasmall=raw.data[,1:3]
  raw.data=raw.data[,4:ncol(raw.data)]
  raw.data=t(raw.data)  
  colnames(raw.data)=metasmall$V1
  }
  # p=stringr::str_split(x,'_',simplify = T)[,2]
  p=stringr::str_split(stringr::str_split(x,'-',simplify = T)[,1],"/",simplify = T)[,2]
  print(p)
  sce <- CreateSeuratObject(counts = raw.data,project = p )
  ct=GetAssayData(object = sce, assay = "RNA", slot = "counts") 
  ct=as.matrix(ct)
  
  dir.create(file.path(stringr::str_split(x,'-',simplify = T)[,1]))
  dir.create(file.path(paste0(stringr::str_split(x,'-',simplify = T)[,1],'/filtered_feature_bc_matrix/')))#filtered_feature_bc_matrix 
  data.table::fwrite(data.frame(rownames(ct),rownames(ct),"Gene Expression"),file = paste0(stringr::str_split(x,'-',simplify = T)[,1],'/filtered_feature_bc_matrix/features.tsv.gz'),
              quote = F,sep = '\t',
              col.names = F,row.names = F)
  #genes.tsv = symbol,symbol,GE, features= ensg,symbol,GE
  data.table::fwrite(data.frame(colnames(ct)),file = paste0(stringr::str_split(x,'-',simplify = T)[,1],'/filtered_feature_bc_matrix/barcodes.tsv.gz'),quote = F,sep = '\t',
              col.names = F,row.names = F)
  file=paste0(stringr::str_split(x,'-',simplify = T)[,1],'/filtered_feature_bc_matrix/matrix.mtx')
  sink(file)
  cat("%%MatrixMarket matrix coordinate integer general\n")
  cat("%\n")
  cat(paste(nrow(ct),ncol(ct),sum(ct>0),"\n")) 
  sink()
  
  #再写入表达量信息
  # tmp1=do.call(rbind,lapply(1:ncol(ct),function(i){
  #   return(data.frame(row=1:nrow(ct),
  #                     col=i,
  #                     exp=ct[,i]))
  # }) )
  tmp=data.table::rbindlist(lapply(1:ncol(ct),function(i){
    return(data.frame(row=1:nrow(ct),
                      col=i,
                      exp=ct[,i]))
  }))
  tmp=tmp[tmp$exp>0,]
  head(tmp)
  data.table::fwrite(tmp,file = paste0(stringr::str_split(x,'-',simplify = T)[,1],'/filtered_feature_bc_matrix/matrix.mtx'),quote = F,sep = '\t',
              col.names = F,row.names = F,append = T)
  R.utils::gzip(file = paste0(stringr::str_split(x,'-',simplify = T)[,1],'/filtered_feature_bc_matrix/matrix.mtx'),remove=TRUE,overwrite=T)
})
# 处理成三文件重质控

# #
# adult_pancreas_2020.rds
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\source data\Adult Pancreas dataset\adult_pancreas_2020.rds)")#
fs=r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\adultpancreas2020\)"
# chronic_pancreatitis.rds
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\source data\Chronic Pancreatitis dataset\chronic_pancreatitis.rds)")#
fs=r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\chronicpancreatitis\)"
# fs
sce$sample.ident=sce$sample_ID
sce.all.list <- SplitObject(sce, split.by = "sample.ident")
sce.all.list

# fs=list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\source data\Adult Pancreas dataset\)",pattern = '\\.rds',full.names = T)
# fs

sceList <-  lapply(1:length(names(sce.all.list)), function(y){
  x=paste0(fs,names(sce.all.list)[[y]])
  ct=GetAssayData(object = sce.all.list[[y]], assay = "RNA", slot = "counts") 
  ct=as.matrix(ct)
  print(x)
  dir.create(file.path(stringr::str_split(x,'\\.',simplify = T)[,1]))
  dir.create(file.path(paste0(stringr::str_split(x,'\\.',simplify = T)[,1],'/filtered_feature_bc_matrix/')))#filtered_feature_bc_matrix 
  data.table::fwrite(data.frame(rownames(ct),rownames(ct),"Gene Expression"),file = paste0(stringr::str_split(x,'\\.',simplify = T)[,1],'/filtered_feature_bc_matrix/features.tsv.gz'),
                     quote = F,sep = '\t',
                     col.names = F,row.names = F)
  #genes.tsv = symbol,symbol,GE, features= ensg,symbol,GE
  data.table::fwrite(data.frame(colnames(ct)),file = paste0(stringr::str_split(x,'\\.',simplify = T)[,1],'/filtered_feature_bc_matrix/barcodes.tsv.gz'),quote = F,sep = '\t',
                     col.names = F,row.names = F)
  file=paste0(stringr::str_split(x,'\\.',simplify = T)[,1],'/filtered_feature_bc_matrix/matrix.mtx')
  sink(file)
  cat("%%MatrixMarket matrix coordinate integer general\n")
  cat("%\n")
  cat(paste(nrow(ct),ncol(ct),sum(ct>0),"\n")) 
  sink()
  
  #再写入表达量信息
  # tmp1=do.call(rbind,lapply(1:ncol(ct),function(i){
  #   return(data.frame(row=1:nrow(ct),
  #                     col=i,
  #                     exp=ct[,i]))
  # }) )
  tmp=data.table::rbindlist(lapply(1:ncol(ct),function(i){
    return(data.frame(row=1:nrow(ct),
                      col=i,
                      exp=ct[,i]))
  }))
  tmp=tmp[tmp$exp>0,]
  head(tmp)
  data.table::fwrite(tmp,file = paste0(stringr::str_split(x,'\\.',simplify = T)[,1],'/filtered_feature_bc_matrix/matrix.mtx'),quote = F,sep = '\t',
                     col.names = F,row.names = F,append = T)
  R.utils::gzip(file = paste0(stringr::str_split(x,'\\.',simplify = T)[,1],'/filtered_feature_bc_matrix/matrix.mtx'),remove=TRUE,overwrite=T)
})




library(scCancer)
# dataPath <- "./data/GSM5032772"     # The path of cell ranger processed data
# savePath <- "./results/GSM5032772"  # A path to save the results
# sampleName <- "GSM5032772"          # The sample name
setwd(r"(H:\download.big.ac.cn\gsa\CRA001160\other_files\)")
dataPath <- "./PDAC/T1"     # The path of cell ranger processed data
savePath <- "./results/T1"  # A path to save the results
sampleName <- "T1"          # The sample name

library(scCancer)
stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  sampleName = sampleName,
  mix.anno = c(human = "hg38", mouse = "mm10"),
  # authorName = authorName
)


dataPath <- "./data/GSM5032772"      # The path of cell ranger processed data
statPath <- "./results/GSM5032772"   # The path of the scStatistics results
savePath <- "./results/GSM5032772"   # A path to save the results
sampleName <- "GSM5032772"           # The sample name
authorName <- "z"            # The author name to mark the report

anno.results <- runScAnnotation(
  dataPath = dataPath,
  statPath = statPath,
  savePath = savePath,
  bool.runMalignancy = F,
  genome = "hg38",
  authorName = authorName,
  sampleName = sampleName,
  geneSet.method = "average"       # or "GSVA"
)
##################################
library(scCancer)

dataPath <- "./data/GSM5032773"     # The path of cell ranger processed data
savePath <- "./results/GSM5032773"  # A path to save the results
sampleName <- "GSM5032773"          # The sample name
authorName <- "z"           # The author name to mark the report

stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  sampleName = sampleName,
  mix.anno = c(human = "hg38", mouse = "mm10"),
  authorName = authorName
)


dataPath <- "./data/GSM5032773"      # The path of cell ranger processed data
statPath <- "./results/GSM5032773"   # The path of the scStatistics results
savePath <- "./results/GSM5032773"   # A path to save the results
sampleName <- "GSM5032773"           # The sample name
authorName <- "z"            # The author name to mark the report

anno.results <- runScAnnotation(
  dataPath = dataPath,
  statPath = statPath,
  savePath = savePath,
  bool.runMalignancy = F,
  genome = "hg38",
  authorName = authorName,
  sampleName = sampleName,
  geneSet.method = "average"       # or "GSVA"
setwd(r"(H:\download.big.ac.cn\gsa\CRA001160\other_files\PDAC)"))


##############

single.savePaths <- c("./results/GSM5032771", "./results/GSM5032772", "./results/GSM5032773")
sampleNames <- c("GSM5032771", "GSM5032772", "GSM5032773")    # The labels for all samples
savePath <- "./results/sample123-comb"       # A path to save the results
combName <- "sample123-comb"                 # A label of the combined samples
authorName <- "z"                # The author name to mark the report
comb.method <- "Harmony"#"NormalMNN"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")

comb.results <- runScCombination(
  single.savePaths = single.savePaths, 
  sampleNames = sampleNames, 
  savePath = savePath, 
  combName = combName,
  bool.runMalignancy = F,
  genome = "hg38",
  authorName = authorName,
  comb.method = comb.method
)

# RenameCells
# add.cell.id