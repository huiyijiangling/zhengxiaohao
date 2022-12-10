load(file=r"(C:\Users\zxh\Desktop\R\meta collect\tcga_pcawg_rnaseq.Rdata)")
load(r"(C:\Users\zxh\Desktop\R\meta collect\ensg_symbol_relation.Rdata)")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# 
# d = read.csv(your_csv_file)
# ## assume 1st column is ID
# ## 2nd column is FC
# 
# ## feature 1: numeric vector
# geneList = d[,2]
# 
# ## feature 2: named vector
# names(geneList) = as.character(d[,1])
# 
# ## feature 3: decreasing orde
# geneList = sort(geneList, decreasing = TRUE)


enrichWP(gene, organism = "Homo sapiens") 
gseWP(geneList, organism = "Homo sapiens")
x <- enrichPathway(gene=de, pvalueCutoff = 0.05, readable=TRUE)
y <- gsePathway(geneList,   pvalueCutoff = 0.2, pAdjustMethod = "BH", verbose = FALSE)
data(gcSample)
str(gcSample) 
compareCluster
ck <- compareCluster(geneCluster = gcSample, fun = "enrichPathway")
#"enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"
#"groupGO",为啥没哟别的啊坐等更新
library(ReactomePA)
data(geneList, package="DOSE"
)
de <- names(geneList)[abs(geneList) > 1.5]
head(de)
mydf <- data.frame(Entrez=de[1:200],
                   logFC = 100:(-99),
                   group = c(rep('A',100),rep('B',100)),
                   othergroup = c(rep('A',100),rep('B',100)))
mydf=dplyr::arrange(mydf,desc(logFC))
xx.formula <- compareCluster(logFC~group, data=mydf,
                             fun='gseGO', OrgDb='org.Hs.eg.db')

# a formula of type Entrez | logFC ~ group for "gseGO", "gseKEGG" and "GSEA".
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 
#serif times new roman arial
theme_big_simple1 <- function(){ 
  theme_bw(base_size=16, base_family="") %+replace% 
    theme(
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      legend.title=element_text(vjust = 1,face="bold",family="sans",size=24),
      legend.text=element_text(vjust = 1,face="bold",family="sans",size=20),
      legend.justification = "right",
      legend.key.width = unit(40, "pt"),
      axis.line = element_line(color = "black", size = 1, linetype = "solid"),
      axis.ticks = element_line(colour = "black", size = 1),
      panel.background  = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),#panel.grid.major = element_line(color="lightgrey"),
      panel.border = element_blank(),
      # legend.position = ("bottom"),
      plot.title = element_text(face="bold",family="sans",size=24, hjust = 0.0, vjust = 1.75),
      axis.text.x = element_text(face="bold",family="sans",hjust = 0.5,color="black", size=20, margin = margin(t = 4, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(face="bold",family="sans",hjust = 0.5,color="black", size=20, margin = margin(t = 0, r = 4, b = 0, l = 0)),
      axis.title.x = element_text(face="bold",family="sans",hjust = 0.5,margin = margin(t = 10, r = 0, b = 0, l = 0), angle = 0, size = 24),
      axis.title.y = element_text(face="bold",family="sans",hjust = 0.5,margin = margin(t = 0, r = 10, b = 0, l = 0), angle = 90, size = 24),
      axis.ticks.length = unit(0.20, "cm"),
      strip.background = element_rect(color="black", size=1, linetype="solid"),
      strip.text.x = element_text(face="bold",family="sans",size = 20, color = "black"),
      strip.text.y = element_text(face="bold",family="sans",size = 20, color = "black")
    )}


dotplot(ck,label_format = function(x) stringr::str_wrap(x, width=70))+coord_fixed(ratio = 0.8)+ theme_big_simple1()#

dotplot(formula_res)

browseKEGG(kk, 'hsa04110')#kk result
library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))


mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(formula_res)
##           Cluster         group othergroup       ID
## 1 downregulated.A downregulated          A hsa04974
## 2 downregulated.A downregulated          A hsa04510
## 3 downregulated.A downregulated          A hsa04151
## 4 downregulated.A downregulated          A hsa04512
## 5 downregulated.B downregulated          B hsa03320
## 6   upregulated.A   upregulated          A hsa04110
##                        Description GeneRatio  BgRatio       pvalue     p.adjust
## 1 Protein digestion and absorption    16/321 103/8142 2.332583e-06 6.531232e-04
## 2                   Focal adhesion    20/321 201/8142 1.226205e-04 1.716687e-02
## 3       PI3K-Akt signaling pathway    28/321 354/8142 3.256458e-04 3.039361e-02
## 4         ECM-receptor interaction    11/321  88/8142 6.382523e-04 4.467766e-02
## 5           PPAR signaling pathway      5/43  75/8142 4.246461e-05 6.709409e-03
## 6                       Cell cycle    20/220 126/8142 1.233602e-10 3.133350e-08
##         qvalue
## 1 6.212037e-04
## 2 1.632788e-02
## 3 2.890820e-02
## 4 4.249417e-02
## 5 6.302643e-03
## 6 2.804822e-08
##                                                                                                                                             geneID
## 1                                                                 1281/50509/1290/477/1294/1360/1289/1292/1296/23428/1359/1300/1287/6505/2006/7373
## 2                                         55742/2317/7058/25759/56034/3693/3480/5159/857/1292/3908/3909/63923/3913/1287/3679/7060/3479/10451/80310
## 3 55970/5618/7058/10161/56034/3693/4254/3480/4908/5159/1292/3908/2690/3909/8817/9223/4915/3551/2791/63923/3913/9863/3667/1287/3679/7060/3479/80310
## 4                                                                                          7058/3693/3339/1292/3908/3909/63923/3913/1287/3679/7060
## 5                                                                                                                         9370/5105/2167/3158/5346
## 6                                                  4171/993/990/5347/701/9700/898/23594/4998/9134/4175/4173/10926/6502/994/699/4609/5111/1869/1029
##   Count
## 1    16
## 2    20
## 3    28
## 4    11
## 5     5
## 6    20
14.3 Visualization of functional profile comparison
14.3.1 Dot plot
We can visualize the result using the dotplot() method.

dotplot(ck)

dotplot(formula_res)

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)
browseKEGG(kk, 'hsa04110')

cnetplot(ck)#good


library(DOSE)
data(geneList)
de = names(geneList)[1:100]
x = enrichDO(de)
y <- mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

library(ggplot2)
library(forcats)
library(enrichplot)

ggplot(y, showCategory = 20, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle("Enriched Disease Ontology")
# Visualizing rich factor of enriched terms using lolliplot.
# Figure 16.1: Visualizing rich factor of enriched terms using lolliplot.
# 
# A very similar concept is Fold Enrichment, which is defined as the ratio of two proportions, (k/n) / (M/N). Using mutate to add the fold enrichment variable is also easy:
#   
#   mutate(x, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
# 16.5 slice
# We can use slice to choose rows by their ordinal position in the enrichment result. Grouped result use the ordinal position with the group.
# 
# In the following example, a GSEA result of Reactome pathway was sorted by the absolute values of NES and the result was grouped by the sign of NES. We then extracted first 5 rows of each groups. The result was displayed in Figure 16.2.
####################################
library(ReactomePA)
x <- gsePathway(geneList)


y <- arrange(x, abs(NES)) %>% 
  group_by(sign(NES)) %>% 
  slice(1:5)

library(forcats)
library(ggplot2)
library(ggstance)
library(enrichplot)

ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalues), showCategory=10) + 
  geom_col(orientation='y') + 
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL)




# compareCluster

geneClusters	
a list of entrez gene id. Alternatively, a formula of type Entrez~group or a formula of type Entrez | logFC ~ group for "gseGO", "gseKEGG" and "GSEA".

fun	
One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" . Users can also supply their own function.

data	
if geneClusters is a formula, the data from which the clusters must be extracted.

source_from	
If using a custom function in "fun", provide the source package as a string here. Otherwise, the function will be obtained from the global environment.

...	
Other arguments.

# font.size=14)
# digits=2, labs.size=5,
# font.size=14, xlab="", ylab="")

# cex_label_category = 1.2) 
# p2 <- cnetplot(edox, node_label="gene", 
#                cex_label_gene = 0.8) 

# p1 <- dotplot(y, label_format = 20) 
# p2 <- dotplot(y, label_format = function(x) stringr::str_wrap(x, width=20))
# cowplot::plot_grid(p1, p2, ncol=2, labels=c("A", "B")) 