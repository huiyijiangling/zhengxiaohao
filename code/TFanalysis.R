
New_loom <- SeuratDisk::as.loom(sce, filename = "./CP.loom", verbose = T)
New_loom


#  https://resources.aertslab.org/cistarget/
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs

#  https://resources.aertslab.org/cistarget/
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")


# https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/
# https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/
dir.create("cisTarget_databases");
rm(list = ls()) 
library(SCENIC)

# ## Load data
# loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
# library(SCopeLoomR)
# loom <- open_loom(loomPath)
# exprMat <- get_dgem(loom)
# cellInfo <- get_cell_annotation(loom)
# close_loom(loom)
# 
# dim(exprMat)
# exprMat[1:4,1:4]
# head(cellInfo)
# table(cellInfo$CellType)
# library(SCENIC)
# # 保证 cisTarget_databases 文件夹下面有下载好2个1G的文件 "hgnc",
# scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases", nCores=16)#mgi
# 
# # scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"


#构建分析数据
exprMat <- as.matrix(sce@assays$RNA@counts)#表达矩阵
exprMat[1:4,1:4]#查看数据
# cellInfo <- sce@meta.data[,c("celltype","nCount_RNA","nFeature_RNA")]
# colnames(cellInfo) <- c('CellType', 'nGene' ,'nUMI')
# head(cellInfo)
# table(cellInfo$CellType)
#构建scenicOptions对象，接下来的SCENIC分析都是基于这个对象的信息生成的
scenicOptions <- initializeScenic(org = "hgnc", dbDir = "cisTarget_databases", nCores = 4)

genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat.filtered <- exprMat[genesKept, ]
exprMat.filtered[1:4,1:4]
runCorrelation(exprMat.filtered, scenicOptions)
exprMat.filtered.log <- log2(exprMat.filtered + 1)
library(doParallel)
runGenie3(exprMat.filtered.log, scenicOptions)

# Build and score the GRN
exprMat_log <- log2(exprMat + 1)
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)#, coexMethod=c("top5perTarget")) # Toy run settings

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions,exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") #

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
