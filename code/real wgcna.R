## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2018-12-20 15:43:52
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-12-20  First version
### Update Log: 2019-09-10 基于R version 3.5.1 (2018-07-02)
### ---------------
rm(list = ls())  ## 魔幻操作，一键清空~ fpkm等于已经标准化完成了效果同TMM voom 这个是已经删除完了
library(WGCNA) #注意跑这个不要乱加载东西
library(stringr)
library(doParallel)
library(DT)
library(plyr)
library(gplots)
enableWGCNAThreads(nThreads = 8)
options(stringsAsFactors = F)
# load("C:/Users/zxh/Desktop/R/second/STAD_mrna_ENSG_fpkm.Rdata")
# load("C:/Users/zxh/Desktop/R/second/mRNA_fpkm_symbol_ebv.Rdata")
#load("rnatpmlog2001_Liver.Rdata") #原始350
load("rnatpmlog2001_Liver_pca.Rdata")
load("./rnatpm_Liver_WGCNA.Rdata")
# load("TcgaTargetGTEX_phenotype_Liver.Rdata") 
# setwd("./nn2/")
# dat=exprSet
# dat=mRNA_fpkm_symbol_ebv
# rownames(rnatpmlog2001_Liver)=substr(rownames(rnatpmlog2001_Liver),1,15)
# dat=rnatpmlog2001_Liver[which(rownames(rnatpmlog2001_Liver) %in% rownames(rnatpm_Liver_mrna)),]#dat=rnatpm_Liver_mrna
dat=rnaExpr[rownames(PC_0_F),]
# 每次都要检测数据
dat[1:4,1:4] 
dim(dat)
# 这里直接删除离群 只有一个离群值时高高挂起，两个就是单独垂下，其实为了好看可以吧plot(,hang=-1) cex=1, 调整间距 ,leaflab="none" 不显示标签
# dat=subset(dat,select=c(-`GTEX-U8XE-2126-SM-5CHRC`,-`GTEX-WFON-1126-SM-5CHTG`))
# TcgaTargetGTEX_phenotype_Liver=TcgaTargetGTEX_phenotype_Liver[-which(rownames(TcgaTargetGTEX_phenotype_Liver) %in% c("GTEX-U8XE-2126-SM-5CHRC","GTEX-WFON-1126-SM-5CHTG")),]
#
group_list=ifelse(TcgaTargetGTEX_phenotype_Liver624$sample_type=="PrimaryTumor",'tumor','normal')
table(group_list) #table函数，查看group_list中的分组个数

# 接下来走WGCNA流程： 
library(WGCNA)
## WGCNA-step 1 : 构建输入数据，需要严格符合格式，否则会陷入调试代码的坑。
if(T){
  fpkm=dat
  fpkm[1:4,1:4]
  datTraits=data.frame(subtype=group_list)
  datTraits$gsm=colnames(dat)
  head(datTraits)
  table(datTraits$subtype)
  RNAseq_voom <- fpkm 
  ## 因为WGCNA针对的是基因进行聚类，而一般我们的聚类是针对样本用hclust即可，所以这个时候需要转置
  WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:1000],])
  datExpr0 <- WGCNA_matrix  ## top 5000 mad genes
  datExpr <- datExpr0 
  ## 下面主要是为了防止临床表型与样本名字对不上
  sampleNames = rownames(datExpr);
  traitRows = match(sampleNames, datTraits$gsm)
  rownames(datTraits) = datTraits[traitRows, 'gsm']
  
}
#检验缺失的基因 无
if(T){
  gsg = goodSamplesGenes(datExpr, verbose = 3);
  gsg$allOK
}
#检验离群样本 outliers 无
if(T){
  sampleTree = hclust(dist(datExpr), method = "average")
  pdf("kk.pdf",width = 80,height = 100)
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  dev.off()
}
#设置cutheight切之
if(F){
  clust = cutreeStatic(sampleTree, cutHeight = 400, minSize = 10)
  clust
  table(clust)
  #clust
  #0  1 
  #2 98 
  keepSamples = (clust==1)
  datExpr = datExpr[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  # save(datExpr, file = "FPKM-01-dataInput.RData")
}
#检验离群样本 outliers 牛逼的检测用z.k
if(T){
  suppressMessages(library(flashClust))
  # sample network based on squared Euclidean distance note that we
  # transpose the data
  A = adjacency(t(datExpr), type = "distance")
  # this calculates the whole network connectivity
  k = as.numeric(apply(A, 2, sum)) - 1
  # standardized connectivity
  Z.k = scale(k)
  # Designate samples as outlying if their Z.k value is below the threshold
  thresholdZ.k = -2  # often -2.5 可以用-0.4，这样删掉11/625个样本，但实际上没有必要，
  
  # the color vector indicates outlyingness (red)
  outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")
  table(outlierColor)
  # calculate the cluster tree using flahsClust or hclust
  suppressMessages(library(flashClust))
  sampleTree = flashClust(as.dist(1 - A), method = "average")
  # Convert traits to a color representation: where red indicates high
  # values
  datTraits_0=as.numeric(as.factor(datTraits$subtype))
  traitColors = data.frame(numbers2colors(datTraits_0, signed = FALSE))#原为datTraits
  dimnames(traitColors)[[2]] = paste(names(datTraits_0), "C", sep = "")
  datColors = data.frame(outlierC = outlierColor, traitColors)
  # Plot the sample dendrogram and the colors underneath.
  pdf("kk2.pdf",width = 80,height = 80)
  plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors, 
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
  # plot of chunk unnamed-chunk-8
  # 
  # Caption: Cluster tree of mouse liver samples. The leaves of the tree correspond to the mice. The first color band underneath the tree indicates which arrays appear to be outlying. The second color band represents body weight (red indicates high value). Similarly, the remaining color-bands color-code the numeric values of physiologic traits.
  # 
  # The Figure shows the resulting cluster tree where each leaf corresponds to a mouse sample. The first color-band underneath the tree indicates which samples appear outlying (colored in red) according to a low value. The mouse labelled F2_221 is highly outlying (Z.k<-5), which is why we remove it from the subsequent analysis. The other color bands color-code physiological traits. Note that the outlying samples do not appear to have systematically different physiologic trait values. Although the outlying samples are suspicious, we will keep them in the subsequent analysis.
  # 
  # Remove outlying samples from expression and trait data
  remove.samples = Z.k < thresholdZ.k | is.na(Z.k)
  table(remove.samples)
  # the following 2 lines differ from what is written in the book
  datExpr = datExpr[!remove.samples, ]
  datTraits=as.data.frame(datTraits_0)
  datTraits = datTraits[!remove.samples, ]
  datTraits=data.frame(subtype=ifelse(datTraits==1,"normal","tumor"))
  datTraits$gsm=rownames(datExpr)
  ## 下面主要是为了防止临床表型与样本名字对不上
  sampleNames = rownames(datExpr);
  traitRows = match(sampleNames, datTraits$gsm)
  rownames(datTraits) = datTraits[traitRows, 'gsm']
}
## WGCNA-step 2 
## 挑选最佳阈值，在  sft$powerEstimate
if(T){
  powers = c(1:10)#c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut=0.8, verbose = 5)
  #设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
  pdf("WGCNA-step2-beta-value.pdf",width = 8,height = 6)
  # Plot the results:
  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.900;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.900,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}
sft$powerEstimate
sft[["fitIndices"]][["mean.k."]]
## WGCNA-step3 构建加权共表达网络（Weight co-expression network)
## 首先是一步法完成网络构建
if(T){
  net = blockwiseModules(
    datExpr,
    power = sft$powerEstimate,
    maxBlockSize = ncol(datExpr),
    TOMType = "unsigned", minModuleSize = 30,
    reassignThreshold = 0, mergeCutHeight = 0.6,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = F,
    saveTOMFileBase = "AS-green-FPKM-TOM",
    verbose = 3
  )
  table(net$colors) 
}

## 然后是分布法完成网络构建，仅供有探索精神的同学挑战。

## 构建加权共表达网络分为两步：
## 1. 计算邻近值，也是就是两个基因在不样品中表达量的表达相关系数(pearson correlation rho)，
## 参考 2.b.2 in https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
## 2. 计算topology overlap similarity (TOM)。 WGCNA认为，只通过计算两个基因的表达相关系数构建共表达网络是不足够的。
## 于是他们用TOM表示两个基因在网络结构上的相似性，即两个基因如果具有相似的邻近基因，这两个基因更倾向于有相互作用。
## 参考 2.b.3 in https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf


if(F){
  #(1)网络构建 Co-expression similarity and adjacency 
  adjacency = adjacency(datExpr, power = sft$powerEstimate) 
  #(2) 邻近矩阵到拓扑矩阵的转换，Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency);
  dissTOM = 1-TOM
  # (3) 聚类拓扑矩阵 Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  # Plot the resulting clustering tree (dendrogram)
  sizeGrWindow(12,9)
  ## 这个时候的geneTree与一步法的 net$dendrograms[[1]] 性质类似，但是还需要进行进一步处理
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04);
  #(4) 聚类分支的修整 dynamicTreeCut 
  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = 30;
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  table(dynamicMods)
  #4. 绘画结果展示
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  #sizeGrWindow(8,6)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  #5. 聚类结果相似模块的融合，Merging of modules whose expression profiles are very similar
  #在聚类树中每一leaf是一个短线，代表一个基因，
  #不同分之间靠的越近表示有高的共表达基因，将共表达极其相似的modules进行融合
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  #sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  #选择有75%相关性的进行融合
  MEDissThres = 0.25
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  
}

## WGCNA-step 4 ，看看模块颜色
if(T){
  
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  table(mergedColors)
  moduleColors=mergedColors
  # Plot the dendrogram and the module colors underneath
  pdf("WGCNA-step4-genes-modules.pdf",width = 8,height = 6)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  ## assign all of the gene to their corresponding module 
  ## hclust for the genes.
}
# 如果感兴趣可以看看样本的形状分布情况。
# 一般来说，normal和tumor肯定是泾渭分明的。
if(F){
  #明确样本数和基因
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  #首先针对样本做个系统聚类
  datExpr_tree<-hclust(dist(datExpr), method = "average")
  par(mar = c(0,5,2,0))
  plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)
  ## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
  #针对前面构造的样品矩阵添加对应颜色
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$subtype)), 
                                  colors = c( "red","green"),signed = FALSE)
  ## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目
  #  sample_colors <- numbers2colors( datTraits ,signed = FALSE)
  ## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，当然，这样给的颜色不明显，意义不大
  #10个样品的系统聚类树及性状热图
  par(mar = c(1,4,3,1),cex=0.8)
  
  pdf("sample-subtype-cluster.pdf",width = 8,height = 6)
  plotDendroAndColors(datExpr_tree, sample_colors,
                      groupLabels = colnames(sample),
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
}
# intramodularConnectivity
## WGCNA-step 5 
## 这一步主要是针对于连续变量，如果是分类变量，需要转换成连续变量方可使用
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  datTraits$subtype=as.factor(datTraits$subtype)
  design=model.matrix(~0+ datTraits$subtype)
  colnames(design)=levels(datTraits$subtype)
  design
  moduleColors <- labels2colors(net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  pdf("WGCNA-step5-Module-trait-relationships.pdf",width = 8,height = 12)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
}
# 因为本例子是分类变量，而且就是泾渭分明的normal和tumor，所以这个时候大家可以思考一下。
# 这个时候的WGCNA和差异分析的区别是什么。

## WGCNA-step 6

if(T){
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  ## 算出每个模块跟基因的皮尔森相关系数矩
  ## MEs是每个模块在每个样本里面的
  ## datExpr是每个基因在每个样本的表达量
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  
  design
  ## 只有连续型性状才能只有计算
  ## 这里把是否是GC变量0,1进行数值化
  GC = as.data.frame(design[,2]);
  names(GC) = "GC"
  geneTraitSignificance = as.data.frame(cor(datExpr, GC, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(GC), sep="");
  names(GSPvalue) = paste("p.GS.", names(GC), sep="");
  
  module = "turquoise"
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  pdf("WGCNA-step6-Module_membership-gene_significance.pdf",width = 8,height = 6)
  #sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for GC",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  # abline(v=0.8,lwd=1,col="black")
  # abline(h=0.2,lwd=1,col="black")
  dev.off()
}

## WGCNA-step 7 ，对TOM矩阵进行可视化
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  geneTree = net$dendrograms[[1]]; 
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); 
  plotTOM = dissTOM^sft$powerEstimate; 
  diag(plotTOM) = NA; 
  #TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  nSelect = 400#400
  # For reproducibility, we set the random seed
  set.seed(10);
  select = sample(nGenes, size = nSelect);
  selectTOM = dissTOM[select, select];
  # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select];
  # Open a graphical window
  sizeGrWindow(9,9)
  # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
  # the color palette; setting the diagonal to NA also improves the clarity of the plot
  plotDiss = selectTOM^sft$powerEstimate;
  diag(plotDiss) = NA;
  
  pdf("WGCNA-step7-Network-heatmap.pdf",width = 8,height = 6)
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes",col=colorpanel(250,'red',"orange",'lemonchiffon'))
  dev.off()
  
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  ## 只有连续型性状才能只有计算
  ## 这里把是否属 GC 表型这个变量0,1进行数值化
  GC = as.data.frame(design[,2]);
  names(GC) = "GC"
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, GC))
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5);
  
  par(cex = 0.9)
  pdf("WGCNA-step7-Eigengene-dendrogram.pdf",width = 8,height = 6)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                        = 90)
  dev.off()
  
  # Plot the dendrogram
  sizeGrWindow(6,6);
  par(cex = 1.0)
  ## 模块的进化树
  pdf("WGCNA-step7-Eigengene-dendrogram-hclust.pdf",width = 8,height = 6)
  plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                        plotHeatmaps = FALSE)
  dev.off()
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  ## 性状与模块热
  
  pdf("WGCNA-step7-Eigengene-adjacency-heatmap.pdf",width = 8,height = 6)
  plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90)
  dev.off()
  
}

## WGCNA-step 8；导出感兴趣模块的基因
if(T){
  # Select module
  module = "turquoise";
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因
  head(moduleColors);head(probes)
  # 这个时候的 moduleColors 和  probes 是一一对应的。
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  head(modProbes)
}

## WGCNA-step 9 ：导出指定的模块内部的基因之间的关联情况。
if(T){
  # Recalculate topological overlap
  TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); 
  # Select module
  module = "turquoise";
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  ## 也是提取指定模块的基因名
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  ## 模块对应的基因关系矩
  cyt = exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
    nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.02,
    nodeNames = modProbes, 
    nodeAttr = moduleColors[inModule]
  );
  # 导出指定的模块内部的基因之间的关联情况。
  head(cyt$edgeData)
  head(cyt$nodeData)
}
## 对模块的基因集批量GO/KEGG富集分析

FilterGenes_spe <- ((abs(geneModuleMembership[moduleGenes, column]) > 0.8)
                    & (abs(geneTraitSignificance[moduleGenes, 1]) > 0.2) )
table(FilterGenes_spe)
trait_hubGenes_spe <- modProbes[FilterGenes_spe]
hubGene=as.data.frame(trait_hubGenes_spe)
if(length(hubGene)!=0){write.csv(trait_hubGenes_spe,file='trait_hubGenes_spe.csv',quote = T)
}else{
  print("no hub under 0.8 MM 0.2 GS")}

