rm(list = ls()) 
gc()
library(Matrix)
# library(scCancer)
library(cowplot);library(dplyr);library(ggExtra);library(ggplot2);library(grid);library(gridExtra);library(GSVA)
library(knitr);library(markdown);library(Matrix);library(methods);library(NNLM);library(pheatmap);library(R.utils)
library(reshape2);library(scds);library(Seurat);library(SingleCellExperiment);library(SoupX)
library(stringr);library(survival);library(survminer);library(harmony);library(liger);library(DropletUtils)
# setwd(r"(C:\Users\zxh\Desktop\R\scCancer)")
source("C:/Users/zxh/Desktop/R/scCancer/scR221/utils.R")
source("C:/Users/zxh/Desktop/R/scCancer/scR221/cnvFunction.R")
source("C:/Users/zxh/Desktop/R/scCancer/scR221/scStatistics.R")
source("C:/Users/zxh/Desktop/R/scCancer/scR221/scAnnotation.R")
source("C:/Users/zxh/Desktop/R/scCancer/scR221/scCombination.R")

options(Matrix.warnDeprecatedCoerce=0)#否则matrix版本报错
source("C:/Users/zxh/Desktop/R/scCancer/scR221/runScAnnotation_unless.R")
source("C:/Users/zxh/Desktop/R/scCancer/scR221/scCombination_unless_preprocess.R")

library(Seurat)
# 
library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 16000 * 1024^2)#4g*32，老师报错


runScAnnotation
CRA001160 t
GSE202051 t
GSE154778 t
GSE156405 t
GSE201333 n
GSE137766 n

GSE134355 n
GSE84133 n

## 系统报错改为英文
Sys.setenv(LANGUAGE="en")
# anno.filter = c("mitochondrial", "ribosome", "dissociation")

# 大概25min
# 查看分析报表，即全部分析结果。根据报表里的解释也能很清楚做了哪些分析，得到哪些结论。下面会有具体解读~
#   
#   （3）其它输入文件情况
# 我们处理自己或公共的数据集有时单细胞数据可能并不会提供完全的两套数据，而是
# 
# 
# 
# 情况1：只有filtered_feature_bc_matrix(~1w cells)
# 由于raw_feature_bc_matrix只起辅助作用，所以没有太大影响；而且也没有解决办法，如果作者没提供的话。
# 
# 情况2：只有raw_feature_bc_matrix(~50w cells)
# 此时我们可以提前加载DropletUtils包，runScStatistics()函数会自动调用该包，过滤raw_feature_bc_matrix里的empty droplet，再进行计算过滤指标操作。
# 
# 情况3：只有单细胞表达矩阵(一般是过滤empty droplet后的count矩阵)信息
# 我们可以使用scCancer包提供的generate10Xdata()函数根据表达矩阵生成filtered_feature_bc_matrix文件夹，然后进行上述的分析。
# 
# 情况4：多个样本的合并
# 首先对每个样本进行上述的两部分析；
# 然后使用scCombination()函数进行合并校正，涉及多种校正算法。示例代码如下
# 
# raw：直接合并，不做校正
# SeuratMNN：调用Seurat包(V3)的MNN合并算法
# NormalMNN：考虑到肿瘤恶性细胞异质性，故结合runScAnnotation()注释结果，仅对非肿瘤恶性细胞进行MNN合并校正【默认校正方法】
# Regression：对假设的系统性的批次效应进行回归校正
# LIGER： 调用liger包的合并算法
# Harmony：调用harmony包的合并算

#integrate metastatistic

GSE158356 #肝5个
GSE156405 #转4个
GSE154778 #肝1个其他1个
meta_pre=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\metaandprecancerlist.xlsx)",col_names = F)


library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)
# stopCluster(cl)
#
# meta pre ipmn
# GSE158356
folder=r"(C:\Users\zxh\Desktop\R\scCancer\GSE158356\)"
fs=list.files(folder,pattern = "GSM")#
fs
x=fs
# GSE156405
folder=r"(C:\Users\zxh\Desktop\R\scCancer\GSE156405\)"
fs=list.files(folder,pattern = "GSM")#
fs
x=fs
# GSE154778
folder=r"(C:\Users\zxh\Desktop\R\scCancer\GSE154778\)"
fs=list.files(folder,pattern = "GSM")#
fs
x=fs
# GSE165399
folder=r"(C:\Users\zxh\Desktop\R\scCancer\GSE165399\)"
fs=list.files(folder,pattern = "GSM")#
fs
x=fs
# PRJEB40416 gastric 
folder=r"(C:\Users\zxh\Desktop\R\scCancer\PRJEB40416\)"
fs=list.files(folder,pattern = "B|F")#mo
fs
x=fs
#
folder=r"(J:\cra\other_files\PDAC\)"
fs=list.files(folder,pattern = "T|N")#mo
fs
x=fs
#
folder=r"(C:\Users\zxh\Desktop\R\scCancer\GSE125588\)"
fs=list.files(folder,pattern = "GSM")#mo
fs
x=fs
#
folder=r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA\)"
fs=list.files(folder,pattern = "")#mo
fs
x=fs
#
folder=r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\)"
fs=list.files(folder,pattern = "T|U")#mo
fs
x=fs
#
folder=r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\chronicpancreatitis\)"
fs=list.files(folder,pattern = "_")#mo
fs
x=fs
#
folder=r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\adultpancreas2020\)"
fs=list.files(folder,pattern = "_")#mo
fs
x=fs

if(F){
scc <- function(x){

  dataPath <- paste0(folder,x)# The path of cell ranger processed data
  savePath <- paste0(folder,"results\\",x)  # A path to save the results
  sampleName <- x        # The sample name
  authorName <- "Xiaobei"           # The author name to mark the report
  stat.results <- runScStatistics(
    dataPath = dataPath,
    savePath = savePath,
    sampleName = sampleName,
    authorName = authorName
  )
  print(x)
  dataPath <- paste0(folder,x)   # The path of cell ranger processed data
  statPath <- paste0(folder,"results\\",x)  # The path of the scStatistics results
  savePath <- paste0(folder,"results1\\",x)   # A path to save the results
  sampleName <- x          # The sample name

  anno.results <- runScAnnotation(
    dataPath = dataPath,
    statPath = statPath,
    savePath = savePath,
    sampleName = sampleName,
    bool.runDiffExpr = F,#这样快
    bool.runMalignancy = F,
    nCell.min = 3,#comb 样本3，对于单个细胞应该是0
    bool.rmContamination = F,
    geneSet.method = "average",       # or "GSVA"
    bool.filter.cell = T,
    bool.filter.gene = T,
    anno.filter = c("mitochondrial", "ribosome", "dissociation"),
    bgPercent.max = 1,

    vars.add.meta = c("mito.percent", "ribo.percent", "diss.percent"),
    vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent"),
    pc.use = 20,
    resolution = 0.3,
    clusterStashName = "default",
    show.features = NULL,
    bool.add.features = T,

    n.markers = 5,
    species = "human",
    genome = "hg38",
    hg.mm.mix = F,
    bool.runDoublet = T,
    doublet.method = "bcds",
    bool.runCellClassify = F,
    ct.templates = NULL,
    coor.names = c("tSNE_1", "tSNE_2"),

    cnv.ref.data = NULL,
    cnv.referAdjMat = NULL,
    cutoff = 0.1,
    p.value.cutoff = 0.5,
    bool.intraTumor = F,
    bool.runCellCycle = F,
    bool.runStemness = F,
    bool.runGeneSets = F,
    geneSets = NULL,

    bool.runExprProgram = F,
    nmf.rank = 50,
    bool.runInteraction = F,
    genReport = T
  )
  
}}



scc_unless<- function(x){
  
  dataPath <- paste0(folder,x)# The path of cell ranger processed data
  savePath <- paste0(folder,"results\\",x)  # A path to save the results
  sampleName <- x        # The sample name
  authorName <- "Xiaobei"           # The author name to mark the report
  stat.results <- runScStatistics(
    dataPath = dataPath,
    savePath = savePath,
    species = "human",#"mouse",
    mix.anno = c(human = "hg38", mouse = "mm10"),#里面是5000和>3%所以不对，最后再改吧
    sampleName = sampleName,
    authorName = authorName
  )
  rm(stat.results)
  print(x)
  dataPath <- paste0(folder,x)   # The path of cell ranger processed data
  statPath <- paste0(folder,"results\\",x)  # The path of the scStatistics results
  savePath <- paste0(folder,"results\\",x)   # A path to save the results
  sampleName <- x          # The sample name
  authorName <- "Xiaobei"           # The author name to mark the report
  anno.results <- runScAnnotation_unless(
    dataPath = dataPath,
    statPath = statPath,
    savePath = savePath,
    sampleName = sampleName,
    bool.runDiffExpr = F,#这样快
    bool.runMalignancy = T,
    nCell.min = 0,#comb 样本3，对于单个细胞应该是0
    bool.rmContamination = F,
    geneSet.method = "average",       # or "GSVA"
    
    bool.filter.cell = T,
    bool.filter.gene = T,
    anno.filter = c("mitochondrial","ribosome"),#,"dissociation"
    bgPercent.max = 1,
    
    vars.add.meta = c("mito.percent", "ribo.percent", "diss.percent"),
    vars.to.regress = c("mito.percent"),#,"nCount_RNA","nFeature_RNA","ribo.percent",#注意这里写啥不重要因为其实啥也没做
    pc.use = 30,
    resolution = 0.3,
    clusterStashName = "default",
    show.features = NULL,
    bool.add.features = T,
    authorName = authorName,
    n.markers = 5,
    species = "human",#"mouse",#,
    genome = "hg19",#"hg19","hg38"
    hg.mm.mix = F,
    bool.runDoublet = T,
    doublet.method = "bcds",
    bool.runCellClassify = F,
    ct.templates = NULL,
    coor.names = c("tSNE_1", "tSNE_2"),
    
    cnv.ref.data = NULL,
    cnv.referAdjMat = NULL,
    cutoff = 0.1,
    p.value.cutoff = 0.5,
    bool.intraTumor = F,
    bool.runCellCycle = F,
    bool.runStemness = F,
    bool.runGeneSets = F,
    geneSets = NULL,
    
    bool.runExprProgram = F,
    nmf.rank = 50,
    bool.runInteraction = F,
    genReport = T
  )
  rm(anno.results)
}

# foreach(x=fs)%do% scc(x)
library(doParallel)
foreach(x=fs)%do% scc_unless(x)
stopCluster(cl)#.combine=rbind
#多合一吗
if(T){
seurat_multi=readRDS(paste0(paste0(folder,"results\\",x), "/expr.RDS"))
seurat_one <- SplitObject(seurat_multi, split.by = "orig.ident")
names=names(seurat_one)
split_into_seurat_one <- function(i){
  dir.create(paste0(folder,"results\\",names[[i]]))
  saveRDS(seurat_one[[i]],paste0(folder,"results\\",names[[i]],"\\expr.RDS"))
  file.copy(paste0(savePath,"\\geneManifest.txt"),paste0(folder,"results\\",names[[i]],"\\geneManifest.txt"))}

foreach(i=1:35)%do% split_into_seurat_one(i)}
#CRA
if(T){
single.savePaths <- list.files(r"(J:\cra\other_files\PDAC\results\)",full.names = T,pattern = "T|N")
sampleNames <- list.files(r"(J:\cra\other_files\PDAC\results\)",pattern = "T|N")    # The labels for all samples
savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\results\CRA001160_comb)"       # A path to save the results
combName <- "CRA001160_comb"                 # A label of the combined samples
authorName <- "Xiaobei"                # The author name to mark the report
comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
}
# GSE202051
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\)",full.names = T,pattern = "T|U")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\)",pattern = "T|U")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\comb)"       # A path to save the results
  combName <- "comb"                 # A label of the combined samples
  authorName <- "Xiaobei"                # The author name to mark the report
  comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
}
# sctcga
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA\results\)",full.names = T,pattern = "")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA\results\)",pattern = "")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA_all\results\comb)"       # A path to save the results
  combName <- "comb"                 # A label of the combined samples
  authorName <- "Xiaobei"                # The author name to mark the report
  comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
}
# sctcga_ecm
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA\results\)",full.names = T,pattern = "")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA\results\)",pattern = "")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA_storma\results\comb)"       # A path to save the results
  combName <- "comb"                 # A label of the combined samples
  authorName <- "Xiaobei"                # The author name to mark the report
  comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
}
# sctcga_immune
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA\results\)",full.names = T,pattern = "")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA\results\)",pattern = "")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA_immune\results\comb)"       # A path to save the results
  combName <- "comb"                 # A label of the combined samples
  authorName <- "Xiaobei"                # The author name to mark the report
  comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
}
# PRJEB40416
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\PRJEB40416\results)",full.names = T,pattern = "EP")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\PRJEB40416\results)",pattern = "EP")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\PRJEB40416\results\PRJEB40416_comb)"       # A path to save the results
  combName <- "PRJEB40416_comb"                 # A label of the combined samples
  authorName <- "Xiaobei"                # The author name to mark the report
  comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
}
#sNUC
if(T){
  single.savePaths <- c(
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\adultpancreas2020\results\)",full.names = T,pattern = "_"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\chronicpancreatitis\results\)",full.names = T,pattern = "_"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\)",full.names = T,pattern = "T|U")
  )
  sampleNames <- c(
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\adultpancreas2020\results\)",pattern = "_"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\chronicpancreatitis\results\)",pattern = "_"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\)",pattern = "T|U")
  )  # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\results\snuc_NCPPAAD)"       # A path to save the results
  combName <- "snuc_NCPPAAD"                 # A label of the combined samples
  authorName <- "Xiaobei"                # The author name to mark the report
  comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
}
# meta pre ipmn
if(T){
  single.savePaths <- c(
    list.files(r"(J:\cra\other_files\PDAC\results\)",full.names = T,pattern = "T|N"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE158356\results\)",full.names = T,pattern = "GSM"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE156405\results\)",full.names = T,pattern = "GSM"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE154778\results\)",full.names = T,pattern = "GSM"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE165399\results\)",full.names = T,pattern = "GSM")
    )
  sampleNames <- c(
    list.files(r"(J:\cra\other_files\PDAC\results\)",pattern = "T|N"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE158356\results\)",pattern = "GSM"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE156405\results\)",pattern = "GSM"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE154778\results\)",pattern = "GSM"),
    list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE165399\results\)",pattern = "GSM")
  )  # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\results\pre_ipmn_primary_meta_comb)"       # A path to save the results
  combName <- "pre_ipmn_primary_meta_comb"                 # A label of the combined samples
  authorName <- "Xiaobei"                # The author name to mark the report
  comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
}
# single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE125588\results\)",full.names = T,pattern = "GSM")
# sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE125588\results\)",pattern = "GSM")    # The labels for all samples
# savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\GSE125588\results\comb)"       # A path to save the results
# combName <- "comb"                 # A label of the combined samples
# authorName <- "Xiaobei"                # The author name to mark the report
# comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")

# # The paths of all sample's "runScAnnotation" results
# single.savePaths <- list.files("H:\\download.big.ac.cn\\gsa\\CRA001160\\results\\",full.names = T,pattern = "T|N")
# sampleNames <- list.files(r"(H:\download.big.ac.cn\gsa\CRA001160\other_files\results\)",pattern = "T|N")    # The labels for all samples
# savePath <- "./results/comb"       # A path to save the results
# combName <- "comb"                 # A label of the combined samples
# authorName <- "Xiaobei"                # The author name to mark the report
# comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")

# single.savePaths <- c("./results/sample1", "./results/sample1", "./results/sample1")
# sampleNames <- c("sample1", "sample2", "sample3")    # The labels for all samples
# savePath <- "./results/sample123-comb"       # A path to save the results
# combName <- "sample123-comb"                 # A label of the combined samples
# authorName <- "Xiaobei"                # The author name to mark the report
# comb.method <- "NormalMNN"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
# 
# comb.results <- runScCombination(
#   single.savePaths = single.savePaths, 
#   sampleNames = sampleNames, 
#   savePath = savePath, 
#   combName = combName,
#   authorName = authorName,
#   comb.method = comb.method
# )
comb.results <- runScCombination_unless(#runScCombination_unless
  single.savePaths = single.savePaths,
  vars.to.regress = c("mito.percent"), #,"G2M.Score","S.Score","nCount_RNA","nFeature_RNA","ribo.percent","nCount_RNA",
  # "G2M.Score","S.Score" 不好用，我并不想出现Ki67为特征的群
  #"nCount_RNA","nFeature_RNA",不要加这个会引起不合适的亚群融合
  sampleNames = sampleNames, 
  savePath = savePath, 
  combName = combName,
  authorName = authorName,
  pc.use = 30,
  resolution = 0.8,#人cra0.3
  bool.runDiffExpr = F,
  bool.runMalignancy = F,
  species ="human",#"mouse",#
  genome = "hg19",#"mm10",#"hg19",hg38
  # coor.names = c("tSNE_1", "tSNE_2"),
  bool.runCellClassify = F,
  bool.runCellCycle = F,
  # bool.runGeneSets = F,
  bool.runExprProgram = F,
  nmf.rank = 50,
  coor.names = c("UMAP_1", "UMAP_2"),
  comb.method = comb.method
)  

# start run subcelltype
# publish
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\Subcelltype\results)",full.names = T,pattern = "T|U")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\Subcelltype\results)",pattern = "T|U")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\Subcelltype\results\comb)"       # A path to save the results
}
# publish
if(F){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\fibro\Subcelltype\results\)",full.names = T,pattern = "T|N")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\fibro\Subcelltype\results\)",pattern = "T|N")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\fibro\Subcelltype\results\comb)"       # A path to save the results
}
# immune
if(F){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\immune\Subcelltype\results\)",full.names = T,pattern = "T|N")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\immune\Subcelltype\results\)",pattern = "T|N")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\immune\Subcelltype\results\comb)"       # A path to save the results
}
# epithelial
if(F){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\epithelial\Subcelltype\results\)",full.names = T,pattern = "T|N")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\epithelial\Subcelltype\results\)",pattern = "T|N")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\epithelial\Subcelltype\results\comb)"       # A path to save the results
}

if(F){
single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\Subcelltype\results\)",full.names = T,pattern = "T|N")
sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\Subcelltype\results\)",pattern = "T|N")    # The labels for all samples
savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\Subcelltype\results\comb)"       # A path to save the results
}
if(F){
single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE125588\results\Subcelltype\results\)",full.names = T,pattern = "GSM")
sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE125588\results\Subcelltype\results\)",pattern = "GSM")    # The labels for all samples
savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\GSE125588\results\Subcelltype\results\comb)"       # A path to save the results
}



# sNUC epi
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\sNUC\epi\Subcelltype\results\)",full.names = T,pattern = "")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\sNUC\epi\Subcelltype\results\)",pattern = "")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\results\snuc_NCPPAAD_split\epi_comb_tmp)"       # A path to save the results
}

# sNUC imm
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\sNUC\imm\Subcelltype\results\)",full.names = T,pattern = "")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\sNUC\imm\Subcelltype\results\)",pattern = "")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\results\snuc_NCPPAAD_split\imm_comb_tmp)"       # A path to save the results
}

# sNUC sto
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\sNUC\sto\Subcelltype\results\)",full.names = T,pattern = "")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\sNUC\sto\Subcelltype\results\)",pattern = "")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\results\snuc_NCPPAAD_split\sto_comb_tmp)"       # A path to save the results
}



# pre_ipmn_primary_meta_comb epi
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\pre_ipmn_primary_meta_comb\epi\Subcelltype\results\)",full.names = T,pattern = "")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\pre_ipmn_primary_meta_comb\epi\Subcelltype\results\)",pattern = "")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\results\pre_ipmn_primary_meta_comb_split\epi_comb_tmp)"       # A path to save the results
}

# pre_ipmn_primary_meta_comb imm
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\pre_ipmn_primary_meta_comb\imm\Subcelltype\results\)",full.names = T,pattern = "")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\pre_ipmn_primary_meta_comb\imm\Subcelltype\results\)",pattern = "")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\results\pre_ipmn_primary_meta_comb_split\imm_comb_tmp)"       # A path to save the results
}

# pre_ipmn_primary_meta_comb sto
if(T){
  single.savePaths <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\pre_ipmn_primary_meta_comb\sto\Subcelltype\results\)",full.names = T,pattern = "")
  sampleNames <- list.files(r"(C:\Users\zxh\Desktop\R\scCancer\pre_ipmn_primary_meta_comb\sto\Subcelltype\results\)",pattern = "")    # The labels for all samples
  savePath <- r"(C:\Users\zxh\Desktop\R\scCancer\results\pre_ipmn_primary_meta_comb_split\sto_comb_tmp)"       # A path to save the results
}



combName <- "comb"                 # A label of the combined samples
authorName <- "Xiaobei"                # The author name to mark the report
comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")


# library(future)
# plan("multisession", workers = 8)
# options(future.globals.maxSize = 16000 * 1024^2)#4g*32
comb.results <- runScCombination_unless(#runScCombination_unless
  single.savePaths = single.savePaths,
  vars.to.regress = c( "mito.percent"), #"nCount_RNA","ribo.percent"
  sampleNames = sampleNames, 
  savePath = savePath, 
  combName = combName,
  authorName = authorName,
  pc.use = 50,
  resolution = 1,#人cra0.3
  bool.runDiffExpr = F,
  bool.runMalignancy = F,
  species = "human",#"mouse",#"human",
  genome = "hg19",#"mm10",#"hg19","hg38",
  # coor.names = c("tSNE_1", "tSNE_2"),
  bool.runCellClassify = F,
  bool.runCellCycle = F,
  bool.runGeneSets = F,
  bool.runExprProgram = F,
  nmf.rank = 50,
  coor.names = c("UMAP_1", "UMAP_2"),
  comb.method = comb.method
) 

# harmony
pc.use = 30,
resolution = 0.8,

# sce.all=sce.all.filt
# #拆分为 个seurat子对象
# sce.all.list <- SplitObject(sce.all, split.by = "orig.ident")

sce.all.int <- RunHarmony(sce.all.filt,c( "orig.ident" ), plot_convergence = T)








# NormalizeData()
# ScaleData()
# JackStraw()
# FindMarkers()
# FindIntegrationAnchors()
# FindClusters() - if clustering over multiple resolutions
# 例如，要运行并行版本，您只需要设置future 并照常调用FindMarkers()功能。

