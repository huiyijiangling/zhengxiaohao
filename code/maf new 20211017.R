#将文件转换为maf格式
require(maftools) 
# var.annovar.maf = annovarToMaf(annovar = "all_annovar2", Center = 'NA', refBuild = 'mm10', tsbCol = 'Tumor_Sample_Barcode', table = 'refGene',sep = "\t")
# 
# write.table(x=var.annovar.maf,file="var_annovar_maf",quote= F,sep="\t",row.names=F)

var_maf = read.maf(maf ="C:/Users/zxh/Desktop/R/paad-tcga-gtex/maf_tcga_20211017/mutect/TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.gz",)
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
download icgc please!
