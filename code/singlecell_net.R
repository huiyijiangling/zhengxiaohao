rm(list=ls())
gc()
options(stringsAsFactors = F)
library(Seurat)
library("DropletUtils")
BiocManager::install(ask = F,"DropletUtils")
fs=list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE165399)",pattern = 'matrix.txt.gz',full.names = T)

sceList <-  lapply(fs, function(x){
  a=data.table::fread(x)
  a=as.data.frame(a)
  a[1:4,1:4]
  raw.data=a[,-1]
  rownames(raw.data)=a[,1]
  p=stringr::str_split(x,'_',simplify = T)[,2]
  sce <- CreateSeuratObject(counts = raw.data,project = p )
  ct=GetAssayData(object = sce, assay = "RNA", slot = "counts") 
  ct=as.matrix(ct)
  dir.create(file.path(stringr::str_split(x,'_',simplify = T)[,1]))
  dir.create(file.path(paste0(stringr::str_split(x,'_',simplify = T)[,1],'/raw_feature_bc_matrix/')))
  data.table::fwrite(data.frame(rownames(ct),rownames(ct),"Gene Expression"),file = paste0(stringr::str_split(x,'_',simplify = T)[,1],'/raw_feature_bc_matrix/features.tsv.gz'),
                     quote = F,sep = '\t',
                     col.names = F,row.names = F)
  #genes.tsv = symbol,symbol,GE, features= ensg,symbol,GE
  data.table::fwrite(data.frame(colnames(ct)),file = paste0(stringr::str_split(x,'_',simplify = T)[,1],'/raw_feature_bc_matrix/barcodes.tsv.gz'),quote = F,sep = '\t',
                     col.names = F,row.names = F)
  file=paste0(stringr::str_split(x,'_',simplify = T)[,1],'/raw_feature_bc_matrix/matrix.mtx')
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
  data.table::fwrite(tmp,file = paste0(stringr::str_split(x,'_',simplify = T)[,1],'/raw_feature_bc_matrix/matrix.mtx'),quote = F,sep = '\t',
                     col.names = F,row.names = F,append = T)
  R.utils::gzip(file = paste0(stringr::str_split(x,'_',simplify = T)[,1],'/raw_feature_bc_matrix/matrix.mtx'),remove=TRUE,overwrite=T)
})

#预处理 SCP1096
#f
SCP1096_features=data.table::fread(r"(H:\singlecell_net\SCP1096\expression\5f431e0b771a5b0dbe1f3c3e\treateddata_scp.genes.csv)",header = F)
SCP1096_features$V2=SCP1096_features$V1
SCP1096_features$V3="Gene Expression"
dir.create(file.path(r"(H:\singlecell_net\SCP1096\sc_input)"))
dir.create(file.path(paste0(r"(H:\singlecell_net\SCP1096\sc_input)",'/filtered_feature_bc_matrix/')))
data.table::fwrite(data.frame(SCP1096_features),file = paste0(r"(H:\singlecell_net\SCP1096\sc_input)",'/filtered_feature_bc_matrix/features.tsv.gz'),quote = F,sep = '\t',
                   col.names = F,row.names = F)
#b
SCP1096_barcodes=data.table::fread(r"(H:\singlecell_net\SCP1096\expression\5f431e0b771a5b0dbe1f3c3e\treateddata_scp.barcodes.csv)",header = F)
data.table::fwrite(data.frame(SCP1096_barcodes),file = paste0(r"(H:\singlecell_net\SCP1096\sc_input)",'/filtered_feature_bc_matrix/barcodes.tsv.gz'),quote = F,sep = '\t',
                   col.names = F,row.names = F)
R.utils::gzip(filename = r"(H:\singlecell_net\SCP1096\expression\5f431e0b771a5b0dbe1f3c3e\gene_sorted-treated_data_scp.mtx)",destname=paste0(r"(H:\singlecell_net\SCP1096\sc_input)",'/filtered_feature_bc_matrix/matrix.mtx.gz'),remove=F,overwrite=T)

if(F){
#预处理 SCP1089
#f
SCP1089_features=data.table::fread(r"(H:\singlecell_net\SCP1089\expression\5f3b67ca771a5b0de1476f7d\naivedata_scp.genes.csv)",header = F)
SCP1089_features$V2=SCP1089_features$V1
SCP1089_features$V3="Gene Expression"
dir.create(file.path(r"(H:\singlecell_net\SCP1089\sc_input)"))
dir.create(file.path(paste0(r"(H:\singlecell_net\SCP1089\sc_input)",'/filtered_feature_bc_matrix/')))
data.table::fwrite(data.frame(SCP1089_features),file = paste0(r"(H:\singlecell_net\SCP1089\sc_input)",'/filtered_feature_bc_matrix/features.tsv.gz'),quote = F,sep = '\t',
                   col.names = F,row.names = F)
#b
SCP1089_barcodes=data.table::fread(r"(H:\singlecell_net\SCP1089\expression\5f3b67ca771a5b0de1476f7d\naivedata_scp.barcodes.csv)",header = F,encoding = "Latin-1")
data.table::fwrite(data.frame(SCP1089_barcodes),file = paste0(r"(H:\singlecell_net\SCP1089\sc_input)",'/filtered_feature_bc_matrix/barcodes.tsv.gz'),quote = F,sep = '\t',
                   col.names = F,row.names = F)
R.utils::gzip(filename = r"(H:\singlecell_net\SCP1089\expression\5f3b67ca771a5b0de1476f7d\gene_sorted-naivedata_scp.mtx)",destname=paste0(r"(H:\singlecell_net\SCP1089\sc_input)",'/filtered_feature_bc_matrix/matrix.mtx.gz'),remove=F,overwrite=T)
}

SCP1089_mtx=data.table::fread(r"(H:\singlecell_net\SCP1089\expression\5f3b67ca771a5b0de1476f7d\gene_sorted-naivedata_scp.mtx)",header = F)
head(SCP1089_mtx)
########################





library(scCancer)

dataPath <- r"(H:\singlecell_net\SCP1096\data\SCP1096)"     # The path of cell ranger processed data
savePath <- r"(H:\singlecell_net\SCP1096\results\SCP1096)"  # A path to save the results
sampleName <- "SCP1096"          # The sample name
authorName <- "z"           # The author name to mark the report

stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  sampleName = sampleName,
  mix.anno = c(human = "hg38", mouse = "mm10"),
  authorName = authorName
)


dataPath <- r"(H:\singlecell_net\SCP1096\data\SCP1096)"      # The path of cell ranger processed data
statPath <- r"(H:\singlecell_net\SCP1096\results\SCP1096)"  # The path of the scStatistics results
savePath <- r"(H:\singlecell_net\SCP1096\results\SCP1096)"  # A path to save the results
sampleName <- "SCP1096"           # The sample name
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
save(anno.results,file = "anno_data.Rdata")

sce=anno.results
Idents(sce)
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE,
                              min.pct = 0.25, 
                              thresh.use = 0.25)
library(Seurat)
# devtools::install_github('satijalab/seurat-data')
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)





single.savePaths <- r"(H:\singlecell_net\SCP1096\data\SCP1096)"
sampleNames <- c("SCP1096")    # The labels for all samples
savePath <- r"(H:\singlecell_net\SCP1096\results\sample123-comb)"  # A path to save the results
combName <- "sample123-comb"                 # A label of the combined samples
authorName <- "z"                # The author name to mark the report
comb.method <- "Raw"#"NormalMNN"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")

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
