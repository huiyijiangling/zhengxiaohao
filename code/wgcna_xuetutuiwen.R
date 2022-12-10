rm(list = ls())
options(stringsAsFactors = F)
# 首先读取第一个文件做测试
a=read.table('exprSet/GSM2078124_13_Neonates_CD11c_microglia.htseq.txt.gz')
# 然后批量读取文件
tmp=do.call(cbind,lapply(list.files('exprSet/'),function(x){
  return(read.table(file.path('exprSet/',x)) )
}))
a=tmp[,seq(2,34,by=2)]
rownames(a)=tmp[,1]
library(stringr)
tmp=str_split(list.files('exprSet/'),'_',simplify = T)
colnames(a)=tmp[,1]
datTraits = data.frame(gsm=tmp[,1],
                       mouse=tmp[,3],
                       genetype=ifelse(grepl('CD',tmp[,4]),'pos','neg')
                       
)
RNAseq_voom <- log(edgeR::cpm(a+1)) 
## 因为WGCNA针对的是基因进行聚类，而一般我们的聚类是针对样本用hclust即可，所以这个时候需要转置
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
datExpr0 <- WGCNA_matrix  ## top 5000 mad genes
datExpr <- datExpr0 
datExpr[1:4,1:4]
save(datExpr,datTraits,file = 'input.Rdata')

#step2_beta_value

rm(list = ls())
library(WGCNA)
load(file = 'input.Rdata')
datExpr[1:4,1:4]
if(T){
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  table(powers)
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  png("step2-beta-value.png",width = 800,height = 600)
  
  ?verbose
  # Plot the results:
  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}

sft$powerEstimate
save(sft,file = "step2_beta_value.Rdata")

#step3_genes_modules
rm(list = ls())
library(WGCNA)
load(file = 'input.Rdata')
load(file = "step2_beta_value.Rdata")


##Weight co-expression network
if(T){
  net = blockwiseModules(
    datExpr,
    power = sft$powerEstimate,
    maxBlockSize = 6000,
    TOMType = "unsigned", minModuleSize = 30,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = F, 
    verbose = 3
  )
  table(net$colors) 
}


##模块可视化
if(T){
  
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  table(mergedColors)
  moduleColors=mergedColors
  # Plot the dendrogram and the module colors underneath
  png("step3-genes-modules.png",width = 800,height = 600)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  ## assign all of the gene to their corresponding module 
  ## hclust for the genes.
}


if(T){
  #明确样本数和基因
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  #首先针对样本做个系统聚类
  datExpr_tree<-hclust(dist(datExpr), method = "average")
  #针对前面构造的样品矩阵添加对应颜色
  sample_colors1 <- numbers2colors(as.numeric(factor(datTraits$mouse)), 
                                   colors = c("green","blue","red" ),signed = FALSE)
  sample_colors2 <- numbers2colors(as.numeric(factor(datTraits$genetype)), 
                                   colors = c('yellow','black'),signed = FALSE)
  ssss=as.matrix(  data.frame(mouse=sample_colors1,
                              gt=sample_colors2))
  par(mar = c(1,4,3,1),cex=0.8)
  
  png("sample-subtype-cluster.png",width = 800,height = 600)
  plotDendroAndColors(datExpr_tree, ssss,
                      groupLabels = colnames(sample),
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
}


save(net,moduleColors,file = "step3_genes_modules.Rdata")

#step4_Module_trait_relationships
rm(list = ls())
library(WGCNA)
load(file = 'input.Rdata')
load(file = "step2_beta_value.Rdata")
load(file = "step3_genes_modules.Rdata")


table(datTraits$mouse)
table(datTraits[,2:3])
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  design1=model.matrix(~0+ datTraits$mouse)
  design2=model.matrix(~0+ datTraits$genetype)
  EAEpos=sign(datTraits[,2]=='EAE' & datTraits[,3]=='pos')
  Neonatespos=sign(datTraits[,2]=='Neonates' & datTraits[,3]=='pos')
  design=cbind(design1,design2,EAEpos,Neonatespos)
  colnames(design)=c('EAE','Neonates','ctrl','neg','pos','EAE_POS','NEO_POS')
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
  png("step4-Module-trait-relationships.png",width = 800,height = 1200,res = 120)
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
  # 除了上面的热图展现形状与基因模块的相关性外
  # 还可以是条形图,但是只能是指定某个形状
  # 或者自己循环一下批量出图。
  if(F){
    png("step4-Module-single_trait-relationships.png",width = 800,height = 1200,res = 120)
    par(mar = c(6, 8.5, 3, 3));
    layout(matrix(c(rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2),rep(6,2),rep(7,2)),byrow = T,nrow = 7))
    for (i in 1:7) {
      y = as.data.frame(design[,i]);
      GS1=as.numeric(cor(y,datExpr, use="p"))
      GeneSignificance=abs(GS1)
      # Next module significance is defined as average gene significance.
      ModuleSignificance=tapply(GeneSignificance,
                                moduleColors, mean, na.rm=T)
      #sizeGrWindow(8,7)
      #par(mfrow = c(1,1))
      # 如果模块太多，下面的展示就不友好
      # 不过，我们可以自定义出图。
      plotModuleSignificance(GeneSignificance,moduleColors)
    }
  }
  save(design,file = "step4_design.Rdata")
}

#step5_Annotation_Module_GO_BP：
rm(list = ls())
library(WGCNA)
load(file = 'input.Rdata')
load(file = "step2_beta_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")


table(moduleColors)
group_g=data.frame(gene=colnames(datExpr),
                   group=moduleColors)
save(group_g,file='step5_Annotation_Module_GO_BP_group_g.Rdata')

rm(list = ls())
load(file='step5_Annotation_Module_GO_BP_group_g.Rdata')
library(clusterProfiler)
# Convert gene ID into entrez genes
tmp <- bitr(group_g$gene, fromType="ENSEMBL", 
            toType="ENTREZID", 
            OrgDb="org.Mm.eg.db")

de_gene_clusters=merge(tmp,group_g,by.x='ENSEMBL',by.y='gene')
table(de_gene_clusters$group)


#list_de_gene_clusters <- split(de_gene_clusters$ENTREZID,de_gene_clusters$group)


# Run full GO enrichment test
formula_res <- compareCluster(
  ENTREZID~group, 
  data=de_gene_clusters, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

# Run GO enrichment test and merge terms 
# that are close to each other to remove result redundancy
lineage1_ego <- simplify(
  formula_res, 
  cutoff=0.5, 
  by="p.adjust", 
  select_fun=min
)

# Plot both analysis results
#pdf('female_compared_GO_term_DE_cluster.pdf',width = 11,height = 6)
#dotplot(formula_res, showCategory=5)
#dev.off()
pdf('Microglia_compared_GO_term_DE_cluster_simplified.pdf',width = 13,height = 8)
dotplot(lineage1_ego, showCategory=5)
dev.off()


# Save results
#write.csv(formula_res@compareClusterResult, 
#          file="Microglia_compared_GO_term_DE_cluster.csv")
write.csv(lineage1_ego@compareClusterResult, 
          file="Microglia_compared_GO_term_DE_cluster_simplified.csv")

## 扩展：enric其它数据库，KEGG, RECTOME, GO(CC,MF) 
