20220508 crispr list
用来做 overlap的
mhighall5_left #  top 5 high 
survOutput005 #  0.05 
load(r"(C:\Users\zxh\Desktop\R\houseessentialdependency\HRT_gene.Rdata)")
load(r"(C:\Users\zxh\Desktop\R\houseessentialdependency\Wang.Rdata)")
load(r"(C:\Users\zxh\Desktop\R\houseessentialdependency\Homo_HK_2013.Rdata)")
over=list(HRT_V1_human_mapped$gene_id,mhighall5_left,rownames(survOutput005),Wang_ensg$ENSEMBL,Homo_HK_2013$ENSEMBL)
names(over)=c("HRT_V1_human","mhighall5_left","survOutput005","Wang","Eisenberg")
library(VennDiagram)
venn.diagram(over,filename = "over.tiff",
             fill = c("dodgerblue", "goldenrod1", "darkorange1","seagreen3"),
             col = "transparent",
             alpha = 0.50,
             label.col = c("orange", "white", "darkorchid4", "white", 
                           "white", "white", "white", "white", "darkblue", "white", 
                           "white", "white", "white", "darkgreen", "white"),
             cex = 1.5,
             fontfamily = "serif",
             fontface = "bold",
             cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
             cat.cex = 1.5,
             cat.pos = 0,
             cat.dist = 0.07,
             cat.fontfamily = "serif",
             # rotation.degree = 270,
             margin = 0.2,imagetype = "tiff")

venn.diagram(over,filename = "5ge.tiff",height = 10000, width = 10000,resolution =1200,cex=1.5,cat.cex = 2,imagetype = "tiff",
             fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),scaled=T,
             alpha = 0.50,	cat.dist = c(0.2, 0.2, 0.15, 0.22,0.2),cat.pos = c(0,-20,180,180,30),margin = 0.02)


# resultTlist_U=subset(resultTlist_U,resultTlist_U$mhighall5_left==T&resultTlist_U$survOutput005==T)

resultTlist_U=get.venn.partitions(over,keep.elements = T)
resultTlist_U=subset(resultTlist_U,resultTlist_U$HRT_V1_human==T)
resultTlist_U=subset(resultTlist_U,resultTlist_U$mhighall5_left==T)
resultTlist_U=subset(resultTlist_U,resultTlist_U$Wang==T)
resultTlist_U=subset(resultTlist_U,resultTlist_U$survOutput005==T)
resultTlist_U$..values..


venn.diagram(resultALLlist,filename = "1Venn.tiff",height = 10000, width = 10000,resolution =1200,cex=1.5,cat.cex = 2,imagetype = "tiff",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,	cat.dist = c(0.2, 0.2, 0.15, 0.22,0.2),cat.pos = c(0,-20,180,180,30),margin = 0.02)

venn.diagram(resultTlist,filename = "2Venn.tiff",height = 10000, width = 10000,resolution =1200,cex=1.5,cat.cex = 2,imagetype = "tiff",
             fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
             alpha = 0.50,	cat.dist = c(0.2, 0.2, 0.15, 0.22,0.2),cat.pos = c(0,-20,180,180,30),margin = 0.02)
venn.diagram(resultNlist,filename = "3Venn.tiff",height = 10000, width = 10000,resolution =1200,cex=1.5,cat.cex = 2,imagetype = "tiff",
             fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
             alpha = 0.50,	cat.dist = c(0.2, 0.2, 0.15, 0.22,0.2),cat.pos = c(0,-20,180,180,30),margin = 0.02)

