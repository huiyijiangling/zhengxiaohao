library(hpgltools)
# meta <- as.data.frame(fission@colData)
exprSet_input_sub=rnaCounts[,which(colnames(rnaCounts) %in% phe$sample)]
phe=phe[phe$sample %in% colnames(exprSet_input_sub),]
exprSet_input_sub=rnaCounts[,which(colnames(rnaCounts) %in% phe$sample)]

fission_data=exprSet_input_sub
write.csv(fission_data,file = "rnaCounts.csv")
fission_data=read.csv(file = "rnaCounts.csv")
rownames(fission_data)=fission_data$X
fission_data=fission_data[,-1]
# duplicated(metaMatrix.RNA$sample)
meta <- as.data.frame(colnames(fission_data))

## Make conditions and batches
# meta$condition <- paste(meta$strain, meta$minute, sep=".")
# meta$batch <- meta$replicate
meta$sample.id <- meta$`colnames(fission_data)`
## Grab the count data
# fission_data <- fission@assays$data$counts
?normalize_expt
## This will make an experiment superclass called 'expt' and it contains
## an ExpressionSet along with any arbitrary additional information one might want to include.
## Along the way it writes a Rdata file which is by default called 'expt.Rdata'
fission_expt <- create_expt(metadata=meta,count_dataframe=fission_data)
fun_norm_2=log2(fission_data+1)
fun_norm <- normalize_expt(fission_expt,norm = "quant",transform="log2")#里面加1了其实就是log2(x+1)
#sm()
fun_norm_1 <- exprs(fun_norm);dim(fun_norm_1)
boxplot(fun_norm_1)
boxplot(rnaCounts)
boxplot(log2(rnaCounts+1))
boxplot(rna)
normalize_expt(fission_expt)
normalize_expt(
  rnaCounts,
  transform = "raw",
  norm = "raw",
  convert = "raw",
  batch = "raw")