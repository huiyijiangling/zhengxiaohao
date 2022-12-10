#我们也是一个非常准备

rm(list = ls())
options(stringsAsFactors = F)
# 
# library(DOSE)
# library(clusterProfiler)
# library(enrichplot)
# 
# 
# data(geneList)
# de <- names(geneList)[1:100]
# x <- enrichDO(de)
# heatplot(x)
# table(de)
# ?enrichDO
# p=as.data.frame(geneList)
# 
# 
# 
# Induced GO DAG graph
# Gene Ontology (GO) is organized as a directed acyclic graph. An insighful way of looking at the results of the analysis is to investigate how the significant GO terms are distributed over the GO graph. The goplot function shows subgraph induced by most significant GO terms.

library(clusterProfiler)
data(geneList, package="DOSE")
as.data.frame(geneList)
de <- names(geneList)[abs(geneList) > 2]
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
library(enrichplot)
barplot(ego, showCategory=20)#找geo
dotplot(ego, showCategory=30)#找geo

go <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="all")
dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
## remove redundent GO terms
ego2 <- simplify(ego)
cnetplot(ego2, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

upsetplot(ego)

heatplot(ego2, foldChange=geneList)

go2 <- enrichKEGG(de, OrgDb = "org.Hs.eg.db", ont="all")#
dotplot(go2, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

#kk <- gseKEGG(geneList, nPerm=10000) 这里非常重要
kk <- gseKEGG(geneList, nPerm=10000)
ridgeplot(kk)
##############################################################################
gseaplot(kk, geneSetID = 1, by = "runningScore", title = kk$Description[1])

gseaplot(kk, geneSetID = 1, by = "preranked", title = kk$Description[1])

gseaplot(kk, geneSetID = 1, title = kk$Description[1])
##############################################################################
gseaplot2(kk, geneSetID = 1:10, title = kk$Description[1],subplots = 1)




