# library(clusterProfiler)
# ## Not run: 
# data(gcSample)
# xx <- compareCluster(gcSample, fun="enrichKEGG",
#                      organism="hsa", pvalueCutoff=0.05)
# as.data.frame(xx)
# # plot(xx, type="dot", caption="KEGG Enrichment Comparison")
# 
# ## formula interface
# mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467',
#                             '100127206', '100128071'),
#                    group = c('A', 'A', 'A', 'B', 'B', 'B'),
#                    othergroup = c('good', 'good', 'bad', 'bad', 'good', 'bad'))
# xx.formula <- compareCluster(Entrez~group, data=mydf,
#                              fun='groupGO', OrgDb='org.Hs.eg.db')
# as.data.frame(xx.formula)
# 
# ## formula interface with more than one grouping variable
# xx.formula.twogroups <- compareCluster(Entrez~group+othergroup, data=mydf,
#                                        fun='groupGO', OrgDb='org.Hs.eg.db')
# as.data.frame(xx.formula.twogroups)
# 
# ## End(Not run)
# 
# gmt <- 'wikipathways-20210310-gmt-Homo_sapiens.gmt'
# wp <- read.gmt.wp(gmt)
# ewp <- GSEA(geneList, TERM2GENE=wp[,c("wpid",
#                                       "gene")], TERM2NAME=wp[,c("wpid", "name")])
# data(DE_GSE8057)
# xx <- compareCluster(Gene~time+treatment,
#                      data=DE_GSE8057, fun = enricher,
#                      TERM2GENE=wp[,c("wpid", "gene")],
#                      TERM2NAME=wp[,c("wpid", "name")])
# 
# 

#############################################
library(dplyr)
library(clusterProfiler)
data(geneList, package="DOSE")
## fold change > 2 as DE genes
de <- names(geneList)[abs(geneList) > 2]
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="BP",readable=TRUE)
## use simplify to remove redundant terms
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust",
                 select_fun=min)
kk <- gseKEGG(geneList, organism = "hsa")
## downloaded from https://wikipathways-data.wmcloud.org/current/gmt/
  gmt <- 'C:/Users/zxh/Desktop/R/wikipathways-20210810-gmt-Homo_sapiens.gmt'
wp <- read.gmt.wp(gmt)
ewp <- GSEA(geneList, TERM2GENE=wp[,c("wpid",
                                      "gene")], TERM2NAME=wp[,c("wpid", "name")])
library(ChIPseeker)
## the file can be downloaded using 'downloadGSMbedFiles("GSM1295076")'
file <- "C:/Users/zxh/Documents/R/win-library/4.1/ChIPseeker/extdata/GEO_sample_data/GSM1295076_CBX6_BF_ChipSeq_mergedReps_peaks.bed.gz"
gr <- readPeakFile(file)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- seq2gene(gr, tssRegion=c(-1000, 1000),
                  flankDistance = 3000, TxDb)
library(clusterProfiler)
g <- bitr(genes, 'ENTREZID', 'SYMBOL', 'org.Hs.eg.db')
## downloaded from
https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X
encode <- read.gmt("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")
enricher(g$SYMBOL, TERM2GENE=encode)
data(DE_GSE8057)
xx <- compareCluster(Gene~time+treatment,
                     data=DE_GSE8057, fun = enricher,
                     TERM2GENE=wp[,c("wpid", "gene")],
                     TERM2NAME=wp[,c("wpid", "name")])
?compareCluster
########################
ego2 <- filter(ego, p.adjust < 0.001, Count > 10)
dim(ego2)
ego3 <- mutate(ego, richFactor = Count / as.numeric
               (sub("/\\d+", "", BgRatio)))
ewp2 <- arrange(ewp, desc(abs(NES))) %>%
  group_by(sign(NES)) %>%
  slice(1:5)


library(ggplot2)
library(forcats)
ggplot(ego3, showCategory = 10,
       aes(richFactor,
           fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2",
           "#7e62a3"),
  trans = "log10",
  guide=guide_colorbar(reverse=TRUE,
                       order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  DOSE::theme_dose(12) +
  xlab("Rich Factor") +
  ylab(NULL) +
  ggtitle("Biological Processes")

ggplot(ewp2, showCategory=10,
       aes(NES, fct_reorder(Description, NES),
           fill=qvalues)) +
  geom_col() +
  scale_fill_gradientn(colours=c("#b3eebe",
                                 "#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE)) +
  DOSE::theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("WikiPathways")
###########################################################
library(enrichplot)
#
xx2 <- pairwise_termsim(xx)
cnetplot(xx2)
#
xx2 <- pairwise_termsim(xx)
emapplot(xx2)
#
xx2 <- pairwise_termsim(xx)
library(ggstar)
dotplot(xx2,split="treatment",showCategory = 5,font.size =6)+
  facet_grid(.~treatment, scale="free",margins = F) #,margins = "treatment"
dotplot(xx2, shape = TRUE)
dotplot(xx2, group = TRUE)
# dotplot(xx2, x = "GeneRatio", group = TRUE, size = "count")
?facet_grid
