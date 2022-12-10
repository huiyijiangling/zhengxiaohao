# library(clusterProfiler)
# data(geneList, package="DOSE")
de <- names(geneList)[abs(geneList) > 2]
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
library(enrichplot)
goplot(ego)
# 1
barplot(ego, showCategory=20)
# 2
dotplot(ego, showCategory=30)
# 3
go <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="all")
dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
## remove redundent GO terms
ego2 <- simplify(ego)
cnetplot(ego2, foldChange=geneList)
# 5 文章
cnetplot(ego2, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
# 6 文章
upsetplot(ego)
# 7 文章
kk <- gseKEGG(geneList, nPerm=10000)
ridgeplot(kk)
#没有顺序