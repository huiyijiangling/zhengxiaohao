##示例
if(!require(factoextra))install.packages("factoextra")
if(!require(FactoMineR))install.packages("FactoMineR")
library(factoextra)
library(FactoMineR)
# load("rnatpmlog2001_Liver_pca.Rdata")
# load("./rnatpm_Liver_WGCNA.Rdata")
# load("TcgaTargetGTEX_phenotype_Liver.Rdata")
rownames(rnatpmlog2001_Liver)=substr(rownames(rnatpmlog2001_Liver),1,15)
#1
myfpkm=rnatpmlog2001_Liver[which(rownames(rnatpmlog2001_Liver) %in% rownames(rnatpm_Liver_mrna)),]
#2
myfpkm=rnaExpr[rownames(PC_0_F),]
#
probesetvar = apply(myfpkm, 1, var)  #表达变化大的基因
ord = order(probesetvar, decreasing=TRUE)[1:500]   #前200个基因，或者更多
pca = prcomp(t(myfpkm[ord,]), scale=TRUE)
ss=summary(pca)
# 样本聚类
pdf(file='pca.pdf',height = 80,width = 80)
fviz_pca_ind(pca, label="all", habillage=as.factor(TcgaTargetGTEX_phenotype_Liver624$sample_type),
             addEllipses=TRUE, ellipse.level=0.999936657516334, palette = "Dark2")
dev.off()
4sigma 0.999936657516334
6sigma 0.999999998026825

myfpkm=mirExpr
#
probesetvar = apply(myfpkm, 1, var)  #表达变化大的基因
ord = order(probesetvar, decreasing=TRUE)[1:500]   #前500个基因，没见到明显异常,否则删除3个
pca = prcomp(t(myfpkm[ord,]), scale=TRUE)
ss=summary(pca)
# 样本聚类
pdf(file='pca.pdf',height = 80,width = 80)
fviz_pca_ind(pca, label="all", habillage=as.factor(metaMatrix.MIR$sample_type),
             addEllipses=TRUE, ellipse.level=0.9973, palette = "Dark2")
dev.off()


#一个函数帮你完成复杂的pca统计分析
pca <- PCA(datdat[,-ncol(datdat)], graph = FALSE)
#

ind.out=apply(pca$x, 2, function(x) which( abs(x - mean(x)) > (4* sd(x)) ) ) %>% Reduce(union, .) %>%
  print()

col <- as.numeric(as.factor(TcgaTargetGTEX_phenotype_Liver624$sample_type)); col[ind.out] <- "purple"
qplot(pca$x[, 1], pca$x[, 2], color = I(col), size = I(2)) + coord_equal()
summary(dis)[1:10,]

# Read more: http://www.sthda.com/english/wiki/ggplot2-colors
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100))

代码部分，筛选基因也可以参照另一篇文章，而不一定是选取200个变化最大的基因，R筛选基因：



#绘图：
plot(pca$x[,1:2],col=as.numeric(as.factor(TcgaTargetGTEX_phenotype_Liver$sample_type)),pch=rownames(TcgaTargetGTEX_phenotype_Liver))
# plot(pca$x[,1:2],col=rep(c(1,2,3,4,1,2,3,4),each=3),pch=rep(c(16,17),each=12))
#或者3D：
library(scatterplot3d)
scatterplot3d(pca$x[,1:3],color=as.numeric(as.factor(TcgaTargetGTEX_phenotype_Liver624$sample_type)),pch=rownames(TcgaTargetGTEX_phenotype_Liver624))
#scatterplot3d(pca$x[,1:3],color=rep(c(1,2,3,4,1,2,3,4),each=3),pch=rep(c(16,17),each=12))

#不建议全放进去
#1
datdat=data.frame(t(rnatpmlog2001_Liver),as.factor(TcgaTargetGTEX_phenotype_Liver$sample_type))
#2
datdat=data.frame(t(rnaExpr),as.factor(TcgaTargetGTEX_phenotype_Liver624$sample_type))
