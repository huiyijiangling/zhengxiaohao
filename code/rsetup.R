# 1/(1+exp(-x))
#SCpubr
#scRNAtoolVis 难看
#CyTOF workflow
#CATALYST
# https://bvieth.github.io/powsimR/articles/powsimR.html
# rna检验
# install.packages("installr")
# library(installr)
# updateR()
install.packages("rvcheck")
# rvcheck::update_all()
# rvcheck::check_r()
library(rvcheck)
getwd()    #查看当前目录work directory，若从桌面打开默认路径为"C:/Users/asus/Documents"
.libPaths()      #查看现在的R包安装路径 如果重复安装可以删除
# "C:/Users/zxh/Documents/R/win-library/4.0" "C:/Program Files/R/R-4.1.2/library" 
# 作者：guguaihezi
# 链接：https://www.jianshu.com/p/1017b57f8d79
# 来源：简书
# 简书著作权归作者所有，任何形式的转载都请联系作者获得授权并注明出处。

# system.time(foreach(i = 1:10, .combine = rbind) %:%
#               +               foreach(j = 1:10, .combine = c) %dopar% mean(rnorm(N, i, j)))
# user  system elapsed 
# 0.09    0.00    2.14 

rm(list = ls()) 
# usethis::create_github_token()
# usethis::edit_r_environ()
# open git
# git config --global user.email 'student_506@126.com'
# git config --global user.name 'huiyijiangling'
library(gitcreds)
gitcreds::gitcreds_get()
gitcreds::gitcreds_set()
# git config --global user.email "student_506@126.com"
#"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"
options()$repos 
options()$BioC_mirror
# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror

# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
install.packages("rlang")
install.packages("vctrs")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(ask = F,version = '3.16',timeout = 1200)
if (!requireNamespace("rjson", quietly = TRUE))
  install.packages("rjson")
install.packages("jsonlite")
install.packages("DescTools")
BiocManager::install(ask = F,"KEGG.db")
BiocManager::install(ask = F,c("GSEABase","GSVA","clusterProfiler" ))
BiocManager::install(ask = F,"clusterProfiler",ask = F,update = T)
BiocManager::install(ask = F,c("BiocParallel" ))
BiocManager::install(ask = F,c("limma" ))
BiocManager::install(ask = F,c("GEOquery","impute" ))
BiocManager::install(ask = F,c("genefu"))
BiocManager::install(ask = F,c("org.Hs.eg.db"))
BiocManager::install(ask = F,c("hgu133plus2.db" ))
BiocManager::install(ask = F,c('DESeq2','edgeR'),
                     ask = F,update = F)
BiocManager::install(ask = F,c('airway'),
                     ask = F,update = F)
BiocManager::install(ask = F,c('colorhcplot'),
                     ask = F,update = F)
if (!requireNamespace("tidyr", quietly = TRUE))
  install.packages("tidyr")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
install.packages("readxl")
install.packages('prodlim')
install.packages('survcomp')

install.packages('lars')
install.packages('tableone')
install.packages('caret')
BiocManager::install(ask = F,'enrichplot')
install.packages('leaps')
install.packages('scales')
install.packages('RSQLite')
install.packages('bit')
install.packages('mime')
install.packages('RCurl')

BiocManager::install(ask = F,'RTCGA')
BiocManager::install(ask = F,'RTCGA.clinical')
BiocManager::install(ask = F,'RTCGA.mRNA')
install.packages(c("officer","flextable"))#ReporteRs has been rewritten. The new package is officer
if(!require("randomForest")) install.packages("randomForest",update = F,ask = F)

install.packages("writexl")

install.packages('WGCNA')
install.packages(c("FactoMineR", "factoextra"))
install.packages(c( "pheatmap","ggpubr"))
install.packages("ggraph")
install.packages("htmltools")
install.packages("ggstatsplot")
install.packages("hrbrthemes")

install.packages("multcompView")
install.packages("ggfortify")
install.packages("qpcR")
install.packages("bbmle")
install.packages("MuMIn")
install.packages("GGally")
install.packages("ggside")
devtools::install_github("YuLab-SMU/nCov2019")
devtools::install_github('stacyderuiter/s245')


install.packages(c("rlist", "pipeR"))
BiocManager::install(ask = F,c("preprocessCore"),ask = F,update = T)
if (!requireNamespace("annotate", quietly = TRUE))
  BiocManager::install(ask = F,"annotate")
BiocManager::install(ask = F,c("AnnotationDbi"))
BiocManager::install(ask = F,c("digest"))
BiocManager::install(ask = F,"TCGAbiolinks")
BiocManager::install(ask = F,"backports")
BiocManager::install(ask = F,c("maftools" ))
BiocManager::install(ask = F,c("deconstructSigs" ))
BiocManager::install(ask = F,c("BSgenome.Hsapiens.UCSC.hg38" ))
BiocManager::install(ask = F,c("BSgenome.Hsapiens.UCSC.hg19" ),ask = F,update = T)
BiocManager::install(ask = F,c("GDCRNATools" ),ask = F,update = T)
BiocManager::install(ask = F,'regexPipes')
if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if(!require("rtracklayer")) BiocManager::install(ask = F,"rtracklayer")
BiocManager::install(ask = F,'Hmisc')


BiocManager::install(ask = F,"VennDiagram")#install.packages("VennDiagram")
BiocManager::install(ask = F,"reshape2")
BiocManager::install(ask = F,"pasilla")

install.packages(c("devtools","pheatmap","stringr",
                   "survival","survminer",
                   "glmnet","timeROC","ROCR"))
BiocManager::install(ask = F,"Circos")

if(!require("KEGGREST")) BiocManager::install(ask = F,"KEGGREST",update = F,ask = F)
if(!require("UpSetR")) BiocManager::install(ask = F,"UpSetR",update = F,ask = F)
if(!require("tidyverse")) install.packages("tidyverse",update = F,ask = F)
install.packages("corrr")
install.packages("survRM2")

install.packages("compareC")


BiocManager::install(ask = F,"sva")
install.packages("ppcor")
install.packages("magrittr")
BiocManager::install(ask = F,"biomaRt")
BiocManager::install(ask = F,"gage")
install.packages("doParallel")
BiocManager::install(ask = F,"recount")
install.packages("pamr")
BiocManager::install(ask = F,"VGAM",update = F,ask = F)
BiocManager::install(ask = F,"pscl",update = F,ask = F)

require(devtools)
devtools::install_github("kokrah/cbcbSEQ")
require(cbcbSEQ)
# vignette('cbcbSEQIntro', package='cbcbSEQ')#部分已经不能用故用下方可替代
BiocManager::install(ask = F,"PROPER")
## oh yeah, cbcbSEQ's combatMod no longer works
## It looks to me like the voomMod function is missing a is.na() check and so
devtools::install_github("elsayed-lab/hpgltools")
#browseVignettes(package = 'hpgltools')
BiocManager::install(ask = F,'fission',update = F,ask = F)
BiocManager::install(ask = F,'AnnoProbe',update = F,ask = F)
if(F){
  #目前安装不了1
  library(devtools)
  devtools::install_github("jmzeng1314/AnnoProbe")
  #改为
  install.packages("C:/Users/zxh/Documents/GitHub/AnnoProbe",repos=NULL,type="source")
}
#目前安装不了2
library(devtools)
devtools::install_github("shijianasdf/BasicBioinformaticsAnalysisFromZhongShan")#?("lucky")
#改为
install.packages("pacman")
install.packages("tuneR")
install.packages("C:/Users/zxh/Documents/GitHub/BasicBioinformaticsAnalysisFromZhongShan",repos=NULL,type="source")
BiocManager::install(ask = F,"AnnotationHub")#用来换坐标的，用法搜索https://www.jianshu.com/p/717ec1420678
library(AnnotationHub)
install.packages("PRROC")
install.packages("RColorBrewer")
install.packages("DMwR")
install.packages("gPCA")
install.packages("qgraph")
install.packages("propagate")
install.packages("lsa")
#source("https://bioconductor.org/biocLite.R")
BiocManager::install(ask = F,"affy")
BiocManager::install(ask = F,"simpleaffy")
install.packages("hdnom")
library(devtools)

BiocManager::install(ask = F,c( 'oligo' ))
BiocManager::install(ask = F,c( 'oligoData' ))
BiocManager::install(ask = F,c( 'pd.hg.u133.plus.2' ))
BiocManager::install(ask = F,c( 'pd.hg.u95av2' ))
BiocManager::install(ask = F,c( 'lumi' ))
install.packages("pec")
BiocManager::install(ask = F,c( 'ArrayExpress' ),update =T) 
install.packages("doMC", repos="http://R-Forge.R-project.org")#以后有空研究一下哪个效率高吧
#The doSMP package was intended to be the Windows alternative to doMC. I believe it was eventually taken off CRAN because of problems building it on newer versions of GCC.
#The doParallel and doSNOW packages are probably the most popular foreach backends available for Windows.
BiocManager::install(ask = F,c('impute'),
                     ask = F,update = F)
BiocManager::install(ask = F,c( 'pd.hg.u219' ),update =T) 
BiocManager::install(ask = F,c( 'miRNAmeConverter' ),update =T) 
BiocManager::install(ask = F,c( 'miRBaseConverter' ),update =T) 
BiocManager::install(ask = F,c( 'tidygraph' ),update =T) 
BiocManager::install(ask = F,c( 'tweenr' ),update =T)
BiocManager::install(ask = F,c( 'reactome.db' ),update =T)

BiocManager::install(ask = F,c( 'maEndToEnd' ),update =T)
BiocManager::install(ask = F,c( 'marray' ),update =T)

install.packages( 'doParallel' )
install.packages( 'FactoInvestigate' )
install.packages( 'VIM' )
install.packages( 'BWStest' )

BiocManager::install(ask = F,c( 'GEOmetadb' ),update =T)
library(GEOmetadb)

## ----------------------------------------------------
if( !file.exists("GEOmetadb.sqlite") ) {
  demo_sqlfile <- getSQLiteFile(destdir = getwd(), destfile = "GEOmetadb.sqlite.gz", type = "demo")
} else {
  demo_sqlfile <- "GEOmetadb.sqlite"
}

library(devtools)
# install.packages("C:/Users/zxh/Documents/GitHub/zstexplorer",repos=NULL,type="source")

BiocManager::install(ask = F,c( 'fda.usc' ),update =T)

BiocManager::install(ask = F,"pROC",update =T)
install.packages('ggDCA')#https://mp.weixin.qq.com/s/dcN1BvmuSO7osWFPPq3pYg
install.packages("e1071")
install.packages( 'My.stepwise' )

# 
# install.packages( 'c060' )
# 
# packages <- c("ggplot2", "dplyr", "reshape2", "RColorBrewer", "scales", "grid","regclass","caret","glmnet",
#               "pROC","AppliedPredictiveModeling","dplyr","ggpubr","ranger","ggplot2","ggfortify","tidyr","survminer","survival", "qdap", "table1", "xtable", "glmnet", "penalizedSVM", "parallel", "randomForestSRC",
#               "ggRandomForests","xtable","peperr","tgp","mlegp","pamr","lattice","c060","limma")
# 
# LoadRpak(packages)
BiocManager::install(ask = F,c( 'Agi4x44PreProcess' ),update =T)
install.packages( 'plotmo' )
Agi4x44PreProcess
install.packages("do")
install.packages("ggloop")
install.packages("randomcoloR")
BiocManager::install(ask = F,c( 'triebeard' ),update =T)
BiocManager::install(ask = F,"GenomicFeatures")
install.packages("stringi")
install.packages("ipflasso")
install.packages("splitstackshape")
install.packages("nomogramEx")

install.packages("riskRegression")

install.packages("foreach")
install.packages("doFuture")
#安装方法：
# install.packages("devtools")

# devtools::install_github("Github-Yilei/ggcor")
# install.packages("C:/Users/zxh/Downloads/ggcor-hou-master.zip",repos=NULL)
# library(ggcor)

install.packages("MatchIt")

install.packages("captioner")
BiocManager::install(ask = F,c( 'ReactomePA' ),update =T)
install.packages("ggridges")

install.packages("writexl")
BiocManager::install(ask = F,c( 'Synapser' ),update =T)
BiocManager::install(ask = F,c( 'GOSemSim' ),update =T)
BiocManager::install(ask = F,c( 'ChIPseeker' ),update =T)
install.packages("ggnewscale")
install.packages("ggstar")
install.packages("tidygraph")
library(tidygraph)
install.packages("grafify")#画图
install.packages("devtools")

install.packages("ROSE")
BiocManager::install(ask = F,c("igraph"),ask = F,update = T)
BiocManager::install(ask = F,c("shiny","markdown","shinyWidgets","plotly","igraph","visNetwork","tibble","shinyjs","yonder"))
install.packages("gridtext")
install.packages("ggcorrplot")
install.packages("corrplot")
BiocManager::install(ask = F,c( 'CATALYST' ),update =T)
install.packages("DiagrammeR")
install.packages("CeRNASeek")
install.packages("easyPubMed")
install.packages("bigmemory")
install.packages("future")
install.packages("future.apply")
install.packages("RCircos")
install.packages("enc")
install.packages("GenVisR")
install.packages("RIdeogram")
install.packages("KEGGREST")
BiocManager::install(ask = F,c( 'KEGGREST' ),update =T)
install.packages("utf8")
install.packages("devtools")
devtools::install_github("Moonerss/CIBERSORT")
install.packages("oncoPredict")
devtools::install_github("maese005/oncoPredict")

install.packages("magick")
BiocManager::install(ask = F,c( 'EnrichedHeatmap' ),update =T)
BiocManager::install(ask = F,c( 'CancerSubtypes' ),update =T)



depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor", "timeROC","pracma")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(ask = F,depen,update = FALSE)
}

if (!requireNamespace("EPIC", quietly = TRUE))
  devtools::install_github("GfellerLab/EPIC", ref="master")
if (!requireNamespace("estimate", quietly = TRUE)){
  rforge <- "http://r-forge.r-project.org"
  install.packages("estimate", repos=rforge, dependencies=TRUE)
}
if (!requireNamespace("IOBR", quietly = TRUE))
  remotes::install_github("IOBR/IOBR",force = T)

# install.packages("devtools")
# system("git config --global --unset http.proxy")
devtools::install_github("DongqiangZeng0808/Blasso")
credentials::set_github_pat(force_new = FALSE, validate = interactive(), verbose = validate)
install.packages("UCSCXenaTools")
BiocManager::install(ask = F,"phenopath")
BiocManager::install(ask = F,"scatterpie")
install.packages('reticulate')
# BiocManager::install(ask = F,"umap")
devtools::install_github("cole-trapnell-lab/L1-graph")


install.packages("tensorflow")
install.packages("keras")
library(tensorflow)
library(keras)
library(reticulate)
reticulate::install_miniconda()
reticulate::conda_create(envname = "py37",python_version = 3.7)
reticulate::use_condaenv("C:/Users/zxh/AppData/Local/r-miniconda/envs/py37")
reticulate::py_install("diopy",pip=T,envname = 'py37')

# Current channels:
#   
#   - https://conda.anaconda.org/conda-forge/win-64
# - https://conda.anaconda.org/conda-forge/noarch
# - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/win-64
# - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/noarch
# - https://repo.anaconda.com/pkgs/main/win-64
# - https://repo.anaconda.com/pkgs/main/noarch
# - https://repo.anaconda.com/pkgs/r/win-64
# - https://repo.anaconda.com/pkgs/r/noarch
# - https://repo.anaconda.com/pkgs/msys2/win-64
# - https://repo.anaconda.com/pkgs/msys2/noarch
# 
# To search for alternate channels that may provide the conda package you're
# looking for, navigate to
# 
#     https://anaconda.org
# 
# and use the search bar at the top of the page.
# 
# 
# Error: one or more Python packages failed to install [error code 1]


# Verify if it was properly selected
# linux内的我使用的是sc
# conda create -n sc python=3.7 （pytorch 是我自己取的名字）
reticulate::py_config()

# To install Keras
# keras::install_keras(method = 'conda', envname = 'py37',)
# tensorflow::install_tensorflow(method = 'conda', envname = 'py37')

# reticulate::py_install('scikit-learn', pip = F) 

# keras::is_keras_available()
# reticulate::py_install('umap-learn')#, pip = T, pip_ignore_installed = T) # Ensure the latest version of UMAP is installed
# reticulate::py_install("louvain", envname = 'py37')

BiocManager::install(ask = F,c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                               'limma', 'S4Vectors', 'SingleCellExperiment',
                               'SummarizedExperiment', 'batchelor', 'Matrix.utils'),update =T)
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages("pastecs")

install.packages("jiebaR")
install.packages("wordcloud")
install.packages("CatEncoders")
install.packages("mltools")
install.packages("g3viz")
install.packages('heatmaply')
BiocManager::install(ask = F,'GenePattern')
install.packages('useful')
install.packages("aplot")
install.packages("patchwork")
install.packages("cowplot")
install.packages("datapasta")
install.packages("DataExplorer")
install.packages("nortest")
install.packages("RNOmni")
install.packages("openxlsx")
install.packages("ggExtra")


install.packages("clustree")
install.packages("harmony")#open版install_github("immunogenomics/harmony")
install.packages("ipflasso")
devtools::install_github("jackwasey/icd")
install.packages("comorbidity")

install.packages("smcfcs")

#rmarkdown
install.packages('rmarkdown')
install.packages('tinytex')
tinytex::install_tinytex()
install.packages("rmdformats")
install.packages("rticles")
install.packages("prettydoc")
install.packages("tufte")

# install.packages("shinythemes")#cerulean

install.packages("dynplot")


install.packages("tsne")
BiocManager::install(ask = F,'bcellViper')
BiocManager::install(ask = F,'viper')

BiocManager::install(ask = F,'scater')
BiocManager::install(ask = F,c("M3Drop"))
BiocManager::install(ask = F,c("scran"))
BiocManager::install(ask = F,c("progeny"))

devtools::install_local("PATH/TO/DIRECTORY/CytoTRACE_0.3.3.tar.gz")
#https://cytotrace.stanford.edu/
TCSeq
BiocManager::install(ask = F,"TCseq")
# pip install cellphonedb
## 查看网页版说明
# browseVignettes("TCseq")
BiocManager::install(ask = F,"pcaMethods")
devtools::install_github("velocyto-team/velocyto.R")
BiocManager::install(ask = F,'Mfuzz')
BiocManager::install(ask = F,'BioNet')
install.packages("SCENIC")


# import scvelo as scv
# 4、Dynamo分析RNA速率【Cell】

# 5、PHATE对单细胞数据降维【Nature】
# http://events.jianshu.io/p/38a9376f5286
# pip install cellphonedb

## Required
BiocManager::install(ask = F,c("AUCell", "RcisTarget"))
BiocManager::install(ask = F,c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(ask = F,c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(ask = F,c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(ask = F,c("doMC", "doRNG"))#装不上
install.packages("doMC", repos="http://R-Forge.R-project.org")#以后有空研究一下哪个效率高吧
# To export/visualize in http://scope.aertslab.org
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC") 

dbFiles
for(featherURL in dbFiles){
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
  #  (1041.7 MB)
  # 
}

install.packages("magrittr")
install.packages("Greg")
install.packages("Gmisc")

#这些方法可能只是可能方法的一个子集。我知道timereg包有一些非常花哨的时间系数处理。我对包的个人体验是有限的，示例中提供的视觉上不那么吸引人的图表让我感到气馁，并且没有适当的插图来解释细节（上次检查 v 1.8.9）。同样，flexsurv应该能够处理比例风险假设。如果你知道任何其他然后请给我发电子邮件

install.packages("BED")
# 不会用啊
# remotes::install_github("patzaw/BED")
#这是用来 https://rdrr.io/github/patzaw/BED/
# https://rdrr.io/github/patzaw/BED/f/vignettes/BED.Rmd
BiocManager::install(ask = F,c("genomeIntervals"))
install.packages("treemap")
install.packages("ade4")
install.packages("adespatial")
install.packages("timeDate")
BiocManager::install(ask = F,c("DEqMS"))
BiocManager::install(ask = F,"lmer")
install.packages("lmer")
install.packages("sampling")
install.packages("infer")
BiocManager::install(ask = F,c("ALDEx2"))
# install.packages("sctransform") # 以前的已经帮忙装了
install.packages("bench")
remotes::install_github("PYangLab/Cepo", dependencies = TRUE,
                        build_vignettes = TRUE)

devtools::install_github('zhanghfd/contamDE')
BiocManager::install(ask = F,"SingleR")
BiocManager::install(ask = F,"GeneOverlap")
install.packages("pdftools")


# file:///C:/Users/zxh/AppData/Local/R/win-library/4.2/contamDE/doc/Intro.html

#escape enrichIt 可以算分


if(F){
  #安装archr必备
  # gnu gsl
  # gtk cairo
  system('
wget https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz
tar -xzvf gsl-latest.tar.gz
cd gsl-2.7.1
./configure
make
sudo make install
sudo apt-get install libgsl0-dev
sudo apt-get  install cairo* libxt*
sudo apt-get install libcairo2-dev
sudo apt-get install libtiff-dev
# gdal gdal-devel
sudo apt install libgeos-dev
# sudo apt-get install build-essential cmake
# wget https://download.osgeo.org/geos/geos-3.9.3.tar.bz2
# tar -xf geos-3.9.3.tar.bz2
# cd geos-3.9.3
# mkdir build
# cd build
# #指定安装位置
# cmake -DCMAKE_INSTALL_PREFIX=/usr/local/geos ..
# #安装到系统默认的位置
# cmake ..
# make
# sudo make install
')
  ld_path <- paste(Sys.getenv("LD_LIBRARY_PATH"), "/usr/local/lib/", sep = ";")
  Sys.setenv(LD_LIBRARY_PATH = ld_path)
  # https://www.jianshu.com/p/de59dd927066 macs3
}
BiocManager::install("rgeos")
install.packages("Seurat")
install.packages("ggdist")
remotes::install_github('satijalab/seurat-wrappers')
install.packages('Cairo')
BiocManager::install(ask = F,"DirichletMultinomial")
BiocManager::install(ask = F,"TFBSTools")
BiocManager::install(ask = F,"motifmatchr")
BiocManager::install(ask = F,"chromVAR")
devtools::install_github("saezlab/dorothea")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
ArchR::installExtraPackages()

devtools::install_github("xmc811/Scillus", ref = "development")
#https://github.com/xmc811/Scillus 美化软件

# for the current development version:
devtools::install_github("HelenaLC/muscat", ref = "devel")

if(!require("estimate"))install.packages("estimate", repos="http://r-forge.r-project.org",dependencies=TRUE)

install.packages("ISOpureR")
install.packages("SimDesign")
install.packages("mRMRe")

BiocManager::install(ask = F,c("MethylSeekR"))

BiocManager::install(ask = F,"DeMixT")
BiocManager::install(ask = F,c("smotefamily"))
install.packages("ggalluvial")
install.packages("Boruta")
install.packages("xgboost")

remotes::install_github("mojaveazure/seurat-disk")
install.packages("scRNAstat")
# 听说你还缺PBMC单细胞数据
devtools::install_github('satijalab/seurat-data')
SeuratData::InstallData("pbmc3k") # 下载不下来
SeuratData::InstallData("ifnb")
BiocManager::install('TENxPBMCData')
install.packages("dslabs")
install.packages("ctv")
ctv::update.views("Distributions")
#gamlss
install.packages("gamlss.add")
devtools::install_github('oganm/homologene')
#propagate
# ::fitDist	
# https://cran.r-project.org/web/views/Distributions.html


# x <-  c(37.50,46.79,48.30,46.04,43.40,39.25,38.49,49.51,40.38,36.98,40.00,38.49,37.74,47.92,44.53,44.91,44.91,40.00,41.51,47.92,36.98,43.40)
# library(fitdistrplus)
# descdist(x, discrete = FALSE)
# 在此处输入图像描述
# 
# 现在您可以尝试拟合不同的分布。例如：
# 
# normal_dist <- fitdist(x, "norm")
# abs 随后检查合身：
# plot(normal_dist)
# scanpy                
# https://github.com/clara-parabricks/rapids-single-cell-examples

# 三大sc de比较 
# devtools::install_github("neurorestore/Libra")
devtools::install_github("huiyijiangling/Libraavoiddup")


#被我删了普通的Libra
BiocManager::install(ask = F,'MAST')
BiocManager::install(ask = F,'ROTS')
devtools::install_github("kdzimm/hierarchicell")
devtools::install_github("Al-Murphy/reanalysis_scRNA_seq_benchmark")
devtools::install_github('JiekaiLab/dior')
BiocManager::install(ask = F,"DropletUtils")
# https://hbctraining.github.io/scRNA-seq_online/lessons/sc_exercises_integ_marker_identification.html
# https://github.com/hbctraining/scRNA-seq_online
devtools::install_github("andymckenzie/DGCA")
# correlatePairs 相关
# Scran correlatePairs
BiocManager::install(ask = F,c("scran"))
devtools::install_github("CBIIT-CGBB/scCorr")
install.packages("scLink")
BiocManager::install(ask = F,"simpleSingleCell")#correlatePairs 来自scran#这包老啊
devtools::install_github("seriph78/COTAN")
devtools::install_github("Bishop-Laboratory/correlationAnalyzeR")
BiocManager::install(ask = F,"fcoex")
BiocManager::install(ask = F,"simpleSingleCell")
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6547541/ 一个方法 2019
# https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/misc.html
install.packages("fmsb")
devtools::install_github("ricardo-bion/ggradar")
install.packages("CMplot")
install.packages("dygraphs")
install.packages("forecast")
# 安装并加载包
install.packages("egg")
install.packages("hexbin")

# library(hexbin)
# 
# library(RColorBrewer)
# 
# # 生成数据
# 
# x <- rnorm(mean=1.5, 5000)
# 
# y <- rnorm(mean=1.6, 5000)
# 
# # 绘制二维散点图
# 
# bin<-hexbin(x, y, xbins=40)
# 
# my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
# 
# plot(bin, main="", colramp=my_colors, legend=F)

# 安装并加载包

install.packages("likert")

# 安装并加载包


install.packages("PerformanceAnalytics")# 
install.packages("qcc")# 

# library(likert)
# 
# # 使用PISA量表数据
# 
# data(pisaitems)
# 
# items28 <- pisaitems[, substr(names(pisaitems), 1, 5) == "ST24Q"]
# 
# # 绘制条形图
# 
# l28 <- likert(items28)
# 
# summary(l28)
# 
# plot(l28)

install.packages("treemap")

# library(treemap)
# 
# # 生成数据
# 
# group=c("group-1","group-2","group-3")
# 
# value=c(13,5,22)
# 
# data=data.frame(group,value)
# 
# # 绘制树图
# 
# treemap(data, index="group", vSize="value", type="index")

# meta 敏感性分析
install.packages("mvmeta")
install.packages("meta")#metainf
# metainf(metaresult, pooled="fixed")
# forest(metainf(metaresult), comb.fixed=TRUE)
install.packages('rmeta')
install.packages("msigdbr")
# install.packages("ggsci")
# install.packages("RColorBrewer")
install.packages("paletteer")
BiocManager::install(ask = F,"Nebulosa")
devtools::install_github("davidsjoberg/ggsankey")
# install.packages("SCpubr")
cran_packages <- c("colortools",
                   "dplyr",
                   "enrichR",
                   "forcats",
                   "ggbeeswarm",
                   "ggExtra",
                   "ggplot2",
                   "ggrepel",
                   "ggplotify",
                   "ggtext",
                   "Matrix",
                   "patchwork",
                   "purrr",
                   "rlang",
                   "scales",
                   "Seurat",
                   "stringr",
                   "svglite",
                   "tidyr",
                   "viridis")

install.packages(cran_packages)
github_packages <- c("saezlab/liana")
remotes::install_github(github_packages)
BiocManager::install(ask = F,"cytofWorkflow")
devtools::install_github("enblacar/SCpubr", ref = "v1.0.4-dev-stable")
# FeatureScatter(object = scRNA, feature1 = "MS4A1", feature2 = "CD3D")
BiocManager::install(ask = F,"PureCN")#CN 在tumor only中找 cnv https://bmccancer.biomedcentral.com/articles/10.1186/s12885-022-09602-4#Sec2
BiocManager::install(ask = F,"exomePeak2")
BiocManager::install("RBGL") #安装依赖包

devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier")
install.packages("palmerpenguins")
install.packages("showtext")
install.packages("ggbeeswarm")
install.packages("rcartocolor")
devtools::install_github("sajuukLyu/ggunchull", type = "source")
devtools::install_github("junjunlab/scRNAtoolVis")
# BiocManager::install(ask = F,"singleCellTK")
devtools::install_github("compbiomed/singleCellTK@devel")#子功能不错,但可视化非常非常丑
# 对于所有用户
# 如果您直接运行sctkPythonInstallConda()（或 sctkPythonInstallVirtualEnv()），或者创建了自己的 singleCellTK 特定环境（或 venv），那么在开始在工具包中执行任何分析之前激活环境非常重要。否则，工具包可能无法正确定位所需的依赖项，并可能导致意外行为或崩溃。
# 
# ### For Python venv user
# selectSCTKVirtualEnvironment("{NAME of your VENV}")
# ### For conda environment user
# selectSCTKConda("{NAME of your VENV}")
# plotMASTThresholdGenes
# exportSCEtoAnnData	Export a SingleCellExperiment R object as Python annData object
# exportSCEtoFlatFile	Export a SingleCellExperiment object to flat text files
# exportSCEToSeurat	Export data in Seurat object
devtools::install_github("HelenaLC/muscat", ref = "master")#子功能不错,可视化差
# aggregateData	Aggregation of single-cell to pseudobulk data
# calcExprFreqs	calcExprFreqs
install.packages("ggThemeAssist")
devtools::install_github("dreamRs/esquisse")
install.packages("wesanderson")
install.packages("ggthemes")#ggplot2
devtools::install_github('Mikata-Project/ggthemr')#ggplot2
devtools::install_github("ricardo-bion/ggtech")#ggplot2
devtools::install_github("JLSteenwyk/ggpubfigs")#ggplot2

install.packages("ggalluvial")
install.packages("ggtern")
install.packages("ggh4x")

BiocManager::install("hpar", version = "devel")
remotes::install_github("TongZhou2017/ttfriends")
BiocManager::install(ask = F,"glmGamPoi",update = F)
BiocManager::install(ask = F,"cBioPortalData",update = F)

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') 

devtools::install_github('kostkalab/scds',ref="master")

install.packages("adabag")

install.packages("gbm")

devtools::install_github("13308204545/Chord") # 还有使用python的方法的还能更强

devtools::install_github('xiebb123456/AutomaticCellTypeIdentification',dependencies = T)

BiocManager::install(ask = F,c('scmap',  'CHETAH',    'clustifyr', 'org.Mm.eg.db',   ),update = F)
devtools::install_github("zwj-tina/scibetR")
BiocManager::install(ask = F,c("monocle"))
BiocManager::install(ask = F,c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))
devtools::install_github("cole-trapnell-lab/garnett")
devtools::install_github("BatadaLab/scID")
# library(devtools)
# library(SingleCellExperiment)
library(M3Drop)
devtools::install_github("bm2-lab/scLearn")
devtools::install_github("powellgenomicslab/scPred")
devtools::install_github("pcahan1/singleCellNet")
devtools::install_github("grisslab/scClassifR")
devtools::install_github("BatadaLab/scID")
devtools::install_github("omicsCore/scTyper")
BiocManager::install(ask = F,"LoomExperiment")
devtools::install_github("cellgeni/sceasy")
BiocManager::install(ask = F,"bioDist")
# 注意查看rJava教程
devtools::install_github('xiebb123456/AutomaticCellTypeIdentification')#, INSTALL_opts=c("--no-multiarch")
BiocManager::install(ask = F,c('randomForestSRC', 'plsRcox', 'superpc', 'CoxBoost'))
BiocManager::install(ask = F,c('survivalsvm', 'BART'))
BiocManager::install(ask = F,c('CoxBoost'))

packages=c("preprocessCore","pROC","ggplot2","gridExtra","grid","ggpubr","scales","stringr","qvalue","vegan","stringr","tidyverse","pROC","survival","readr","gdata","DescTools","moments","fastcluster","Rfast","dynamicTreeCut","data.table","plyr","dplyr")
sapply(packages,require,character=TRUE)

cran.packages <- c("msigdbr", "dplyr", "purrr", "stringr","magrittr",
                   "RobustRankAggreg", "tibble", "reshape2", "ggsci",
                   "tidyr", "aplot", "ggfun", "ggplotify", "ggridges",
                   "gghalves", "Seurat", "SeuratObject", "methods",
                   "devtools", "BiocManager","data.table","doParallel",
                   "doRNG")

  install.packages(cran.packages, ask = F, update = F)
# install packages from Bioconductor
bioconductor.packages <- c("GSEABase", "AUCell", "SummarizedExperiment",
                           "singscore", "GSVA", "ComplexHeatmap", "ggtree",
                           "Nebulosa")

  BiocManager::install(bioconductor.packages, ask = F, update = F)

if (!requireNamespace("UCell", quietly = TRUE)) {
  devtools::install_github("carmonalab/UCell")
}
if (!requireNamespace("irGSEA", quietly = TRUE)) {
  devtools::install_github("chuiqin/irGSEA")
}
  
devtools::install_github('smorabit/hdWGCNA', ref='dev')
install.packages('pagoda2')
devtools::install_github("Nanostring-Biostats/GeomxTools", 
                         build_vignettes = TRUE, ref = "dev")
devtools::install_github("dtm2451/dittoSeq")
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("zhanghao-njmu/SCP")
if(F){
  py <- Sys.which("python3")
  reticulate:::python_version(py)
  SCP::PrepareVirtualEnv(python = py, pipy_mirror = "https://pypi.tuna.tsinghua.edu.cn/simple", remove_old = TRUE)
  reticulate::virtualenv_python("SCP")
}

devtools::install_github("omarwagih/ggseqlogo")
BiocManager::install("gggenes",ask = F, update = F)
devtools::install_github("TheHumphreysLab/plot1cell")#HaojiaWu/plot1cell
bioc.packages <- c("biomaRt","GenomeInfoDb","EnsDb.Hsapiens.v86","GEOquery","simplifyEnrichment","ComplexHeatmap")
BiocManager::install(bioc.packages,ask = F, update = F)

# dev.packages <- c("chris-mcginnis-ucsf/DoubletFinder","Novartis/hdf5r","mojaveazure/loomR")
# devtools::install_github(dev.packages)
install.packages("My.stepwise", ask = F, update = F)
#
## 下载依赖包
list.of.packages <- c("SingleCellExperiment", "Biobase", "fastglm", "ggplot2", "monocle", "plyr", "RColorBrewer", "ggrepel", "ggridges", "gridExtra", "devtools", "mixtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages,ask = F, update = F)

## 从 GitHub 上下载 GeneSwitches
devtools::install_github("SGDDNB/GeneSwitches")


devtools::install_github("YosefLab/VISION")
devtools::install_github("wu-yc/scMetabolism")

BiocManager::install("scran",ask = F, update = F)
install.packages("iCellR")

# https://zhuanlan.zhihu.com/p/111382806
Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.1")
install.packages("rjags")
devtools::install_github("broadinstitute/infercnv")
