#将文件转换为maf格式
require(maftools)
library(GDCRNATools)
ECM_2020=read.csv(r"(C:\Users\zxh\Desktop\x学习笔记\CAFplusECM\matrisomeproject\MatrisomeDB\human\916f45ca.csv)")
# var.annovar.maf = annovarToMaf(annovar = "all_annovar2", Center = 'NA', refBuild = 'mm10', tsbCol = 'Tumor_Sample_Barcode', table = 'refGene',sep = "\t")
# 
# write.table(x=var.annovar.maf,file="var_annovar_maf",quote= F,sep="\t",row.names=F)
var_maf = read.maf(maf ="C:/Users/zxh/Desktop/R/paad-tcga-gtex/maf_tcga_20211017/mutect/TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.gz")
#概览maf文件
pdf('plotmafSummary_paad.pdf',width = 9.6,height=6)
plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median')
dev.off()
pdf('oncoplot_paad.pdf',width = 9,height=6)
oncoplot(maf = var_maf, top = 30,showTumorSampleBarcodes = F )# fontSize = 2 ,
dev.off()
#绘制箱线图
pdf('laml.titv_paad.pdf',width = 9,height=6)
laml.titv = titv(maf = var_maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
dev.off()
#分析相互关系图
pdf('somaticInteractions_paad.pdf',width = 9,height=6)
somaticInteractions(maf = var_maf, top = 25, pvalue = c(0.05,0.01),showSum = F)
dev.off()
# download icgc please!




laml.plus.gistic2=data.table::fread("C:/Users/zxh/Desktop/R/paad-tcga-gtex/maf_tcga_20211017/mutect/TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.gz")
laml.plus.gistic2$samplename=substr(laml.plus.gistic2$Tumor_Sample_Barcode,1,15)
TCGA_mut_analysis=unique(laml.plus.gistic2$samplename)
table(laml.plus.gistic2$Variant_Classification)
laml.plus.gistic2=subset(laml.plus.gistic2,laml.plus.gistic2$Variant_Classification!="Silent")
#去掉所有silent
laml.plus.gistic2_krasmut=subset(laml.plus.gistic2,laml.plus.gistic2$Hugo_Symbol=="KRAS")
table(laml.plus.gistic2_krasmut$Variant_Classification)
laml.plus.gistic2$KRAS_OR_WT=ifelse(laml.plus.gistic2$samplename %in% laml.plus.gistic2_krasmut$samplename,"KRASmut","WT")

laml.plus.gistic2_krasmut_miss=subset(laml.plus.gistic2_krasmut,laml.plus.gistic2_krasmut$Variant_Classification=="Missense_Mutation")

table(laml.plus.gistic2_krasmut_miss$Variant_Classification)

laml.plus.gistic2_krasmut_miss_onlyone=tidyr::pivot_wider(
  laml.plus.gistic2_krasmut_miss,
  id_cols = "samplename",
  names_from = c("Hugo_Symbol"),
  names_prefix = "",
  # names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "check_unique",
  values_from = "HGVSp_Short",
  values_fill = NA,
  values_fn = function(x) paste0(x,collapse = ",")
)

table(laml.plus.gistic2_krasmut_miss_onlyone$KRAS)
# openxlsx::write.xlsx(laml.plus.gistic2_krasmut,file = "laml.plus.gistic2_krasmut.xlsx")

laml.plus.gistic2_krasmut_miss_onlyone=subset(laml.plus.gistic2_krasmut_miss_onlyone,KRAS %in%c("p.G12D","p.G12V"))#,"p.G12V"

#WT


laml2=subset(laml.plus.gistic2,select = c("samplename","KRAS_OR_WT"))

laml2=laml.plus.gistic2_krasmut_miss_onlyone
laml2=unique(laml2)
rnaCounts_havemutdata=subset(rnaCounts,select=Reduce(intersect,list(colnames(rnaCounts),laml2$samplename)))
rnaCounts_havemutdata=as.data.frame(rnaCounts_havemutdata)
laml2=laml2[which(laml2$samplename%in%Reduce(intersect,list(colnames(rnaCounts),laml2$samplename))),]
colnames(rnaCounts_havemutdata)==laml2$samplename
laml2=laml2[match(colnames(rnaCounts_havemutdata),laml2$samplename),]
colnames(rnaCounts_havemutdata)==laml2$samplename

DEGAll <- gdcDEAnalysis(counts     = rnaCounts_havemutdata, 
                        group      = laml2$KRAS, #KRAS_OR_WT
                        comparison = 'p.G12D-p.G12V', 
                        method     = 'DESeq2',#edgeR
                        n.cores = 4,
                        filter =T)#默认是true，不筛选则选择false
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all', fc = 1.000, pval = 0.05)#
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 1.000, pval = 0.05)#
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 1, pval = 0.05)#其实放到0.05也可以


laml.plus.gistic2_krasmut$HGVSp_Short

laml.plus.gistic2_krasmut$all_effects
deg.all=dePC
class(DEGAll$logFC)
gdcVolcanoPlot_modify <- function (deg.all, fc = 2, pval = 0.01) 
{
  library(ggplot2)
  library(GDCRNATools)
  # library(ggrepel)
  geneList <- deg.all
  geneList$threshold <- c()
  geneList$threshold[geneList$logFC > log(fc, 2) & geneList$FDR < 
                       pval] <- 1
  geneList$threshold[geneList$logFC >= -log(fc, 2) & geneList$logFC <= 
                       log(fc, 2) | geneList$FDR >= pval] <- 2
  geneList$threshold[geneList$logFC < -log(fc, 2) & geneList$FDR < 
                       pval] <- 3
  geneList$threshold <- as.factor(geneList$threshold)
  lim <- max(max(geneList$logFC), abs(min(geneList$logFC))) + 
    0.5
  volcano <- ggplot(data = geneList,label = `symbol`,#new
                    aes(x = `logFC`,y = -log10(`FDR`)))
  volcano + geom_point(aes(color = `threshold`), alpha = 1, 
                       size = 0.8) + xlab("log2(Fold Change)") + ylab("-log10(FDR)") + 
    scale_colour_manual(breaks = geneList$threshold, values = c("1"="red","2"="black","3"="green3")) + xlim(c(-lim, lim)) + geom_vline(xintercept = c(-log(fc, 
                                                                                                                                          2), log(fc, 2)), color = "darkgreen", linetype = 3) + 
    geom_hline(yintercept = -log(pval, 10), color = "darkgreen", 
               linetype = 3) + theme_bw() + theme(axis.line = element_line(colour = "black"), 
                                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.border = element_rect(colour = "black"), panel.background = element_blank()) + 
    theme(legend.position = "none") + theme(axis.text = element_text(size = 14), 
                                            axis.title = element_text(size = 16))+
    ggrepel::geom_text_repel(
      # data = subset(geneList, geneList$FDR < 0.000001 & abs(geneList$logFC) >= 3),
      data = subset(geneList, geneList$symbol%in% ECM_2020$Gene &geneList$FDR < 0.05 & abs(geneList$logFC) >= 2),
      aes(label = `symbol`),
      size = 3,
      box.padding = unit(0.5, "lines"),
      max.overlaps = getOption("ggrepel.max.overlaps", default = 1000),
      point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
}
pdf('gdcVolcanoPlot_paad.pdf',width = 10,height=5)
gdcVolcanoPlot_modify(dePC, fc = 2, pval = 0.05)
gdcVolcanoPlot(dePC, fc = 3, pval = 0.05)
dev.off()

