rm(list=ls())
options(stringsAsFactors = F)
library(TCGAbiolinks)
# GDCquery(project, data.category, data.type, workflow.type,
#          legacy = FALSE, access, platform, file.type, barcode,
#          experimental.strategy, sample.type)
#先从数据库里找到符合各项参数要求的数据
query <- GDCquery(project = "TCGA-ESCA", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  legacy = FALSE,
                  access = "open",
                  # platform =, 
                  experimental.strategy = "RNA-Seq")
#查询
results<-getResults(query)
dim(results)
results[1:5,1:5]
colnames(results)

#再使用命令GDCdownload(）下载
GDCdownload(query,directory = "./mrna_count/",method = "api",files.per.chunk = 5)
data <- GDCprepare(query,directory = "./mrna_count/")
save(data,file = "ESCA_count_querydata.Rdata")
count_data=assay(data)
count_clinical=colData(data)
save(count_data,count_clinical,file = "ESCA_count.Rdata")

#别混起来


rm(list=ls())
options(stringsAsFactors = F)
library(GDCRNATools)
# 2.1 Normalization of HTSeq-Counts data
# ### load RNA counts data
# data(rnaCounts)
# rnaCounts[1:5,1:5]
# 
# ### load miRNAs counts data
# data(mirCounts)
# mirCounts[1:5,1:5]
# 
# ### Normalization of RNAseq data
# rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)
# rnaExpr[1:5,1:5]
# 
# ### Normalization of miRNAs data
# mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)
# mirExpr[1:5,1:5]
# 2.2 Parse and filter RNAseq metadata
# metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
#                                    data.type  = 'RNAseq', 
#                                    write.meta = FALSE)
# metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
# metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
# metaMatrix.RNA[1:5,]
# 2.3 ceRNAs network analysis
# ### Identification of differentially expressed genes ###
# 
# DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
#                         group      = metaMatrix.RNA$sample_type, 
#                         comparison = 'PrimaryTumor-SolidTissueNormal', 
#                         method     = 'limma')
# DEGAll[1:5,]
# 
# ### All DEGs
# deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')
# deALL[1:5,]
# 
# ### DE long-noncoding genes
# deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')
# deLNC[1:5,]
# 
# ### DE protein coding genes
# dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')
# dePC[1:5,]
# 
# ### ceRNAs network analysis of DEGs
# 
# ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
#                           pc          = rownames(dePC), 
#                           lnc.targets = 'starBase', 
#                           pc.targets  = 'starBase', 
#                           rna.expr    = rnaExpr, 
#                           mir.expr    = mirExpr)
# ceOutput[1:5,]
# 
# 
# ### Export ceRNAs network to Cytoscape
# 
# ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 
#                       & ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]
# 
# ###### Export edges
# 
# edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
# edges[1:5,]
# 
# ##### Export nodes
# nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
# nodes[1:5,]

############################################################
3. Case study: TCGA-ESCA
3.1 Download data
# set up directories for downloaded data
project <- 'TCGA-ESCA'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')

### Download RNAseq data
gdcRNADownload(project.id     = 'TCGA-ESCA', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method = 'gdc-client', ## use gdc-client tool to download data
               directory      = rnadir)

### Download miRNAs data
gdcRNADownload(project.id     = 'TCGA-ESCA', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method = 'gdc-client', ## use gdc-client tool to download data
               directory      = mirdir)
3.2 Data organization
### Parse RNAseq metadata
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-ESCA',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

# Filter duplicated samples in RNAseq metadata
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
# Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
metaMatrix.RNA.ESCA=metaMatrix.RNA
save(metaMatrix.RNA.ESCA,file = "metaMatrix.RNA.ESCA.Rdata")
### Parse miRNAs metadata
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-ESCA',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

# Filter duplicated samples in miRNAs metadata
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)
# Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)

### Merge raw counts data
# Merge RNAseq data
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, 
                         organized = FALSE, ## if target data are in folders
                         data.type = 'RNAseq')

# Merge miRNAs data
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir,
                         organized = FALSE, ## if target data are in folders
                         data.type = 'miRNAs')

### TMM normalization and voom transformation
# Normalization of RNAseq data
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = F)

# Normalization of miRNAs data
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)
save(metaMatrix.MIR,metaMatrix.RNA,mirCounts,mirExpr,rnaCounts,rnaExpr,file="ESCA_GDCRNATOOLS.Rdata")
#subtype
# load(file="ESCA_GDCRNATOOLS.Rdata")
a=read.table(file="./case_lists/ebv55.txt",sep='\t',header = F,stringsAsFactors = F)

count_data=t(metaMatrix.MIR)

colnames(count_data)=count_data[4,]
vars <- as.vector(unlist(unique(a)))
count_data <- as.data.frame(count_data) %>%
  dplyr::select(.,one_of(vars))
count_data=t(count_data)

metaMatrix.MIR=as.data.frame(count_data)

count_data=t(metaMatrix.RNA)
colnames(count_data)=count_data[4,]
vars <- as.vector(unlist(unique(a)))
count_data <- as.data.frame(count_data) %>%
  dplyr::select(.,one_of(vars))
count_data=t(count_data)
rownames(count_data)=count_data[,5]
metaMatrix.RNA=as.data.frame(count_data)

count_data=mirExpr
colnames(count_data)=substr(colnames(count_data),1,15)
vars <- as.vector(unlist(unique(a)))
count_data <- as.data.frame(count_data) %>%
  dplyr::select(.,one_of(vars))
mirExpr=as.data.frame(count_data)

count_data=mirCounts

colnames(count_data)=substr(colnames(count_data),1,15)
vars <- as.vector(unlist(unique(a)))
count_data <- as.data.frame(count_data) %>%
  dplyr::select(.,one_of(vars))
mirCounts=as.data.frame(count_data)

count_data=rnaCounts

colnames(count_data)=substr(colnames(count_data),1,15)
vars <- as.vector(unlist(unique(a)))
count_data <- as.data.frame(count_data) %>%
  dplyr::select(.,one_of(vars))
rnaCounts=as.data.frame(count_data)

count_data=rnaExpr

colnames(count_data)=substr(colnames(count_data),1,15)
vars <- as.vector(unlist(unique(a)))
count_data <- as.data.frame(count_data) %>%
  dplyr::select(.,one_of(vars))
rnaExpr=as.data.frame(count_data)

save(metaMatrix.MIR,metaMatrix.RNA,mirCounts,mirExpr,rnaCounts,rnaExpr,file="ESCA_GDCRNATOOLS_ebv.Rdata")

rm(list=ls())
options(stringsAsFactors = F)
load(file="ESCA_GDCRNATOOLS_ebv.Rdata")
library(GDCRNATools)
### Differential gene expression analysis
DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma',
                        filter =F)#默认是true，不筛选则选择false
#data(DEGAll)

# All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all', fc = 0, pval = 1)#

# DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 1)#

# DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 1)#
Volcano plot and Heatmap
#Volcano plot
gdcVolcanoPlot(DEGAll)
# Barplot
gdcBarPlot(deg = deALL, angle = 45, data.type = 'RNAseq')

#Heatmap
#Heatmap is generated based on the heatmap.2() function in gplots package.
degName = rownames(deALL)
gdcHeatmap(deg.id = degName, metadata = metaMatrix.RNA, rna.expr = rnaExpr)#有问题啊
image.png
image.png
3.3 Competing endogenous RNAs network analysis
（ceRNAs network analysis）

### The 3 steps of ceRNAs network analysis:
# Hypergeometric test
# Pearson correlation analysis
# Regulation pattern analysis

### All of the 3 steps can be performed in a single function


### ceRNAs network analysis using internal databases
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)#超几何分布及p值

# load(file="miRNA_2list.Rdata")
# r=mirExpr[miRNA_2list$V2,]
# 不要删除mmirna会影响作图的，可以删除样本，明白吗
# ceOutput3 <- gdcCEAnalysis(lnc         = rownames(deLNC),
#                            pc          = rownames(dePC),
#                            lnc.targets = 'starBase',
#                            pc.targets  = 'starBase',
#                            rna.expr    = rnaExpr,
#                            mir.expr    = r)

### ceRNAs network analysis using user-provided datasets
# load miRNA-lncRNA interactions
# data(lncTarget)
# lncTarget[1:3]
# 
# load miRNA-mRNA interactions 内置的数据库/或者自己的数据库
# data(pcTarget)
# pcTarget[1:3]
# 
# ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
#                           pc          = rownames(dePC), 
#                           lnc.targets = lncTarget, 
#                           pc.targets  = pcTarget, 
#                           rna.expr    = rnaExpr, 
#                           mir.expr    = mirExpr)
# 'mrna_inter_target'
### Network visulization in Cytoscape
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = lnc_inter_target,
                          pc.targets  = mrna_inter_target,
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)#超几何分布及p值
# Filter potential ceRNA interactions
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 & 
                        ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]
save(ceOutput2,ceOutput,file="LNCinterGENEinter.Rdata")
save(ceOutput2,ceOutput,file="LNCinterGENEstar.Rdata")
save(ceOutput2,ceOutput,file="LNCstarGENEinter.Rdata")
save(ceOutput2,ceOutput,file="LNCstarGENEstar.Rdata")
save(ceOutput2,ceOutput,file="LNCinterGENEinter7015001.Rdata")
save(ceOutput2,ceOutput,file="LNCinterGENEinter6015001.Rdata")
save(ceOutput2,ceOutput,file="LNCinterGENEinter5015001.Rdata")
load(file="LNCinterGENEinter7015001.Rdata")
load(file="LNCinterGENEinter.Rdata")
load(file="LNCinterGENEstar.Rdata")
load(file="LNCstarGENEinter.Rdata")
load(file="LNCstarGENEstar.Rdata")
load(file="LNCinterGENEinter6015001.Rdata")
ceOutput2$lncRNAs
ppppp=ceOutput2[ceOutput2$Genes %in% wgcna_mrna2$nodeName,]
qqqqq=ceOutput4[ceOutput4$lncRNAs %in% wgcna_lnc$nodeName,]
table(qqqqq$Genes)
survOutput2 <- gdcSurvivalAnalysis(gene     = x, 
                                   method   = 'coxph', 
                                   rna.expr = rnaExpr, 
                                   metadata = metaMatrix.RNA)
# KM analysis
survOutput1 <- gdcSurvivalAnalysis(gene     = ceOutput2$Genes, 
                                   method   = 'KM', 
                                   rna.expr = rnaExpr, 
                                   metadata = metaMatrix.RNA, 
                                   sep      = 'median')
x=rownames(survOutput1[survOutput1$pValue<0.05,])
# Edges and nodes can be simply imported into Cytoscape 
# for network visualization
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')

write.csv(edges, file='edges.csv',quote = F) ### Network of Cytoscape
write.csv(nodes, file='nodes.csv',quote = F) ### Table of Cytoscape
### Correlation plot on a local webpage
shinyCorPlot(gene1    = rownames(deLNC), 
             gene2    = rownames(dePC), 
             rna.expr = rnaExpr, 
             metadata = metaMatrix.RNA)
image.png
3.4 Other downstream analyses
Univariate survival analysis

# CoxPH analysis
survOutput <- gdcSurvivalAnalysis(gene     = XXXXXXXXXXXXXXXXx, 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA)
# KM analysis
survOutput <- gdcSurvivalAnalysis(gene     = rownames(dePC), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA, 
                                  sep      = 'median')
XXXXXXXXXXXXXXXXx=rownames(survOutput[survOutput$pValue<0.05,])
x2=dePC[XXXXXXXXXXXXXXXXx,]

survOutput
# KM plot on a local webpage by shinyKMPlot
shinyKMPlot(gene = rownames(deALL), rna.expr = rnaExpr, 
            metadata = metaMatrix.RNA)
image.png
3.5 Functional enrichment analysis
All the functional enrichment analyses can be performed in a single function, including:
  
  Gene Ontology (BP, CC, MF) analysis
KEGG pathway analysis
Disease Ontology analysis
The speed was too slow and taked the top 100.
# Gene Ontology (BP, CC, MF) analysis #The speed is too slow and take the top 100.
enrichOutput <- gdcEnrichAnalysis(gene = rownames(deALL)[1:100], simplify = TRUE)
### This step may take a few minutes ###
# Step 1/5: BP analysis done!
# Step 2/5: CC analysis done!
# Step 3/5: MF analysis done!
# Step 4/5: KEGG analysis done!
# Step 5/5: DO analysis done!

#data(enrichOutput)

# Barplot
gdcEnrichPlot(enrichOutput, type = 'bar', category = 'GO', num.terms = 10)
write.csv(enrichOutput, "enrichOutput.csv")

# Bubble plot
gdcEnrichPlot(enrichOutput, type='bubble', category='GO', num.terms = 10)

# KEGG pathway analysis
gdcEnrichPlot(enrichOutput, type = "bar", category = "KEGG", num.terms = 10, bar.color = "dodgerblue")
#bar.color = "chocolate1"

# Disease Ontology analysis
gdcEnrichPlot(enrichOutput, category='DO',type = 'bubble', num.terms = 20)

# View pathway maps on a local webpage
library(pathview)

deg <- deALL$logFC
names(deg) <- rownames(deALL)
pathways <- as.character(enrichOutput$Terms[enrichOutput$Category=='KEGG'])

shinyPathview(deg, pathways = pathways, directory = 'pathview')
