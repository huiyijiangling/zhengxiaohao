# rnatpm_Pancreas  
library(AnnoProbe)
dataset=rnaCounts_Pancreas[,dat$sample]
probe2gene=data.frame(probe_id=rownames(dataset),symbol=rownames(dataset))
source("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/updateName.R")
forthetree=ensg2SymbolUsingNCBImodify(alias=probe2gene$symbol,gene.info.file = "Homo_sapiens.gene_info.gz")
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
if(F){
  forthetree=updateName(probe2gene$symbol)
  table(forthetree$gene_biotype)
  probe2gene=merge(forthetree,probe2gene,by.x="ALIAS",by.y="symbol",all.y=T)
}
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
probe2gene=unique(probe2gene[,c("probe_id","Symbol")])
# probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
colnames(probe2gene)=c("probe_id","symbol")
probes_expr=dataset
genes_expr <- filterEM(probes_expr,probe2gene)
tide=genes_expr
tide=t(scale(t(log2(tide+1)),center=T,scale = F))
psych::describe(tide)
write.table(tide,gzfile("newfile.txt.gz"),sep = "\t",quote = F)

#非常蛋疼的是我们并没有办法找到一个完全符合他们的example
# food=psych::describe(aaaaa)
# psych::describe(qwer)
# aaaaa=read.table("C:/Users/zxh/Downloads/Prat.self_subtract",header=T)
# psych::describe(paad_tide)
# psych::describe(hulu)

tide_phe=read.csv(file="PAAD tide.csv")
tide_phe=merge(tide_phe,merrrr,by.x="Patient",by.y="sample")
tide_phe=as.data.frame(tide_phe)
ui_in_tide_after=subset(ui,select = tide_phe$Patient)
