#将文件转换为maf格式

# var.annovar.maf = annovarToMaf(annovar = "all_annovar2", Center = 'NA', refBuild = 'mm10', tsbCol = 'Tumor_Sample_Barcode', table = 'refGene',sep = "\t")
# 
# write.table(x=var.annovar.maf,file="var_annovar_maf",quote= F,sep="\t",row.names=F)
library(maf)
var_maf = read.maf(maf ="var_annovar_maf")

#概览maf文件

plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median')

oncoplot(maf = var_maf, top = 30, fontSize = 12 ,showTumorSampleBarcodes = F )

#绘制箱线图

laml.titv = titv(maf = var_maf, plot = FALSE, useSyn = TRUE)

plotTiTv(res = laml.titv)

#分析相互关系图

somaticInteractions(maf = var_maf, top = 15, pvalue = c(0.05, 0.1))