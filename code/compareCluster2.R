library(clusterProfiler)
library(ggplot2)
data(gcSample)
str(gcSample) 
ck <- compareCluster(geneCluster = gcSample, fun = enrichKEGG)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="gseKEGG")

head(formula_res)
dotplot(ck)
dotplot(formula_res)
dotplot(formula_res, x="group") + facet_grid(~othergroup)
#理论上套上去就可以了