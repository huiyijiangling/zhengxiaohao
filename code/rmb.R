
m6a_genelist=readxl::read_excel(r"(C:\Users\zxh\Desktop\杜永星老师 提交\m6a.xlsx)")

markers <- list()
markers$gene <- m6a_genelist$gene
markers$E <- m6a_genelist[m6a_genelist$anno=="Erasers",]$gene
markers$W <- m6a_genelist[m6a_genelist$anno=="Writers",]$gene
markers$R <- m6a_genelist[m6a_genelist$anno=="Readers",]$gene
markers$pasteEW <- m6a_genelist$pasteEW[!is.na(m6a_genelist$pasteEW)]
# markers$pasteR <- m6a_genelist$pasteR[!is.na(m6a_genelist$pasteR)]
DimPlot(sce,label = T)
library(UCell)
library(irGSEA)
library(tidyverse)
sce <- irGSEA::irGSEA.score(object = sce, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 4,
                             min.cells = 3, min.feature = 0,
                            geneset = markers,
                            custom = T,
                            msigdb = F,
                            species = "Homo sapiens", category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell","UCell", "singscore", "ssgsea"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')
save(sce,file="r.Rdata")
table(sce$nCount_AUCell)
halfvlnplot <- irGSEA.halfvlnplot(object = sce,
                                  method = "ssgsea",cluster.levels = 0:18,
                                  show.geneset = "E")
                                    # "HALLMARK-INFLAMMATORY-RESPONSE")
halfvlnplot
ridgeplot <- irGSEA.ridgeplot(object = pbmc3k.final,
                              method = "UCell",
                              show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
ridgeplot
# densityheatmap <- irGSEA.densityheatmap(object = pbmc3k.final,
#                                         method = "UCell",
#                                         show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
# densityheatmap
#
scatterplot[[1]] <- irGSEA::irGSEA.density.scatterplot(object = sce,
                                          method = "AUCell",
                                          show.geneset = "E",
                                          reduction = "umap")
scatterplot[[2]] <- irGSEA::irGSEA.density.scatterplot(object = sce,
                                                  method = "AUCell",
                                                  show.geneset = "W",
                                                  reduction = "umap")
scatterplot[[3]] <- irGSEA::irGSEA.density.scatterplot(object = sce,
                                                  method = "AUCell",
                                                  show.geneset = "R",
                                                  reduction = "umap")
plot_grid(plotlist=list(scatterplot[[1]],scatterplot[[2]],scatterplot[[3]]),ncol = 3)
irGSEA.heatmap.plot2 <- irGSEA.heatmap(object = sce, method = "UCell",
                                       show.geneset = NULL)
esi=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\meta collect\esi.xlsx)")
meta_new_rownames=rownames(sce@meta.data)
meta_new=dplyr::full_join(
  sce@meta.data,
  esi,
  by = c("sample.ident"="type"),
  
)
rownames(meta_new)=meta_new_rownames
sce@meta.data=meta_new

#2000 overlap ecm500
color=c("#be0f2d", "#cf273c", "#cb364a", "#eb4035", "#d55f6f", "#f5abb7", "#3d6756", "#4f7764", "#509a80", "#a7d4c3", "#badde9", "#bbd69d", "#8d4b45", "#92523c", "#a05d46", "#dba880", "#f0c986", "#f1e0dc", "#355386", "#3e68a0", "#5091c0", "#94d2ef", "#bcddea", "#cfe3ef", "#2d1c4d", "#594675", "#684e94", "#b186a3", "#d5c0cf", "#e7ddb8", "#4d584c", "#3e563b", "#46776d", "#6a855b", "#b0ce95", "#dfe7d7", "#543c4f", "#634459", "#785059", "#9a7b81", "#ceacb3", "#e3c7c5", "#414445", "#4a615c", "#517177", "#7e8d7a", "#a2ac9e", "#cbcec1", "#45302f", "#4d3434", "#6e5655", "#59473d", "#af9a8b", "#dccfc0")
color=color[1:length(levels(sce$orig.ident))]
names(color)=levels(sce$orig.ident)

ppp=list()
ppp[[1]]=SCpubr::do_DimPlot(sce,group.by = "orig.ident",reduction = "umap",label=T,
                      colors.use = color,repel = T,
                      legend.position = "left")

pdf("umap tcga.pdf",height = 10,width = 10)
print(ppp[[1]])
dev.off()
ppp[[2]]=SCpubr::do_DimPlot(sce,group.by = "hist2",reduction = "umap",label=T,
                            # colors.use = color,
                            repel = T,
                            legend.position = "left")

ppp[[2]]
pdf("umap tcga esi.pdf",height = 10,width = 10)
print(ppp[[2]])
dev.off()
dittoSeq::dittoDimPlot(sce,"sample.ident",reduction = "umap",do.ellipse = T,do.label = T)

dittoSeq::dittoDimPlot(sce,"hist2",reduction = "umap",do.ellipse = T,do.label = T)

cells.use <- colnames(sce[, sce$sample.ident=="PAAD"])
p <- SCpubr::do_DimPlot(sce, 
                        # split.by = "seurat_clusters", 
                        # ncol = 3, 
                        # idents.keep = c("0", "1", "7"),
                        reduction = "umap",
                        colors.use = "Red",
                        # legend.position = "none",
                        cells.highlight = cells.use,
                        font.size = 24)
pdf("Select PAAD.pdf",height = 28,width = 20)
print(p)
dev.off()

print(ppp[[1]])
# CD274 PDCD1
p5=SCpubr::do_FeaturePlot(sce, features = c(markers$gene),enforce_symmetry=F,size,viridis_color_map = "B",ncol = 5)
pdf("isogene .pdf",height = 28,width = 20)
print(p5)
dev.off()



p6=SCpubr::do_FeaturePlot(sce, features = c("ALKBH5"),enforce_symmetry=F,viridis_color_map = "B")
print(p5)
#为了好看, = F
# "PDCD1","VCAN","PTPRC","ACTA2","FAP","FTO","COL1A1","EPCAM",'PECAM1','MME'
pdf("cell cycle umap tcga ecm.pdf",height = 10,width = 10)
print(p5)
dev.off()
pdf("cell cycle umap tcga ecm density.pdf",height = 5,width = 5)
#BD best
p5=SCpubr::do_NebulosaPlot(sce, features = c("CLDN18"),viridis_color_map = "B")#为了好看,enforce_symmetry = F
# "VCAN","PTPRC","ACTA2","FAP","FTO","COL1A1"  ,,'PECAM1','MME'
print(p5)
dev.off()


plot_grid(plotlist = list(ppp[[1]],p,p5,p6),ncol = 2)

pdf("cell cycle umap tcga.pdf",height = 10,width = 10)
print(DimPlot(sce,
              reduction = "umap",
              group.by= "sample.ident",
              pt.size = 1,label = T, repel = TRUE
              # split.by = "Phase"
)+ggsci::scale_colour_igv())
dev.off()


# +ggplot2::geom_mark_ellipse(aes(fill = sample.ident,
#                                 label = sample.ident),
#                             expand = unit(0.5,"mm"),
#                             label.buffer = unit(-5, 'mm'))