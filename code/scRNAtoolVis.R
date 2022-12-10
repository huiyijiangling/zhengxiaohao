library(scRNAtoolVis)

# httest <- system.file("extdata", "htdata.RDS", package = "scRNAtoolVis")
# pbmc <- readRDS(httest)
pbmc=sce
# load markergene
# markergene <- system.file("extdata", "top5pbmc.markers.csv", package = "scRNAtoolVis")
# markers <- read.table(markergene, sep = ',', header = TRUE)
ddd=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\ori\comb\DE-Genes.RDS)")
markers=ddd
library("scales") 
library(ggsci)

# mycol <- hue_pal()(9)
# mycol1 <- pal_npg()(9)

mycol <- hue_pal()(max(pbmc$comb.cluster))
mycol1 <- pal_npg()(max(pbmc$comb.cluster))

AverageHeatmap(object = pbmc,
               markerGene = markers$gene,
               annoCol = TRUE,
               clusterAnnoName = F,
               htRange = c(-2,0,2),
               htCol=c("Blue", "#210C4AFF", "yellow"),#viridis_pal(option = "B")(3),#c( "#000000","#000000","#fdd49e", "red"),
               myanCol = ggsci::pal_igv("default")(max(pbmc$comb.cluster)),
               width = 16,height = 16) +
AverageHeatmap(object = pbmc,
                 markerGene = markers$gene,
                 annoCol = TRUE,
                 myanCol = mycol1)
# c("Blue", "#210C4AFF", "yellow") 适合以0为主的
# c("Blue", "Black", "yellow"), scale后的适合以0为主的
# viridis_pal(option = "B")(3)#从零开始适合大多数是0的
viridis_pal(option = "B")(9)

#######################################2

library(scRNAtoolVis)

httest <- system.file("extdata", "htdata.RDS", package = "scRNAtoolVis")
pbmc <- readRDS(httest)

# add groups
pbmc$groups <- rep(c('stim','control'),each = 1319)
# add celltype
pbmc$celltype <- Seurat::Idents(pbmc)
# load markergene
data("top3pbmc.markers")
# check
head(top3pbmc.markers,3)

  
jjDotPlot(object = pbmc,
            gene = top3pbmc.markers$gene,
            id = 'celltype',
            # anno = F,
            # markerGene = top3pbmc.markers,#finallmarkers结果直接load进行,而且只能选一个，注意互斥
            xtree = F,
            ytree = F,
            rescale = T,#重整为0-1
            rescale.min = 0,#rescale.min = -2,
            rescale.max = 1,#rescale.max = 2,
            midpoint =0.5,#midpoint =0,
            dot.col	= c("white","#990000"),#c('blue','white','red')
            point.geom = F,
            point.shape = 22,
            tile.geom = T,
            tree.pos = 'left',
            same.pos.label = T,
            # yPosition = 10.3,
            split.by = 'groups',
            split.by.aesGroup = T,#T颜色不用于分组Set split.by.aesGroup = T to turn off the colors grouped by groups:
            gene.order = top3pbmc.markers$gene,#rev(top3pbmc.markers$gene),
            # cluster.order = 8:0,
            plot.margin = c(1,1,1,1))# 上下左右距离
################################################################################### 3
  
  # test <- system.file("extdata", "pbmc.markers.csv", package = "scRNAtoolVis")
  # markers <- read.csv(test)

  # plot 这个非常适合有比例的，而直接的稚嫩用单纯的了
  markerVocalno(markers = markers,
                topn = 5,
                labelCol = ggsci::pal_igv("default")(max(pbmc$comb.cluster)))
                # ggsci::pal_aaas()(9))
