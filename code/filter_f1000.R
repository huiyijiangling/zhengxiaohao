if(F){
library(Biobase)
palmieri_medians <- Biobase::rowMedians(as.matrix(genes_expr_GSE32688_rna))
hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")
man_threshold <- quantile(palmieri_medians,seq(0, 1, 0.1),names = F)[5]
median(palmieri_medians)
hist_res <- hist(palmieri_medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)
no_of_samples <-  table(phenoDat_GSE32688_rna$source_name_ch1)
no_of_samples
samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(as.matrix(genes_expr_GSE32688_rna), 1,function(x){sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
idx_man_threshold
palmieri_manfiltered <- subset(genes_expr_GSE32688_rna, idx_man_threshold)
}
# 完整版
# 一下为简化版
filterEx  <- function(datasets,percentage,group){
palmieri_medians <- Biobase::rowMedians(as.matrix(datasets))
percentage=(percentage*100)+1
man_threshold <- quantile(palmieri_medians,seq(0, 1, 0.01),names = F)[percentage]
no_of_samples <-  table(group)
samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(as.matrix(datasets), 1,function(x){sum(x > man_threshold) >= samples_cutoff})
palmieri_manfiltered <- subset(datasets, idx_man_threshold)
print(quantile(palmieri_medians,seq(0, 1, 0.05),names = T))
hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")
plot(hist_res)
abline(v = man_threshold, col = "coral4", lwd = 2)
return(idx_man_threshold)
}
# flist=filterEx(genes_expr_GSE32688_rna,0.4,phenoDat_GSE32688_rna$source_name_ch1)



gdcVoomNormalization_modify <- function(counts) {
  library(edgeR)
  expr = DGEList(counts = counts)
  ## filter out low expression genes
  keepALL <- rowSums(cpm(expr) > 1) >= 0.5*ncol(counts)
  nGenes <- as.numeric(summary(keepALL)[2]) + 
    as.numeric(summary(keepALL)[3])
  nKeep <- summary(keepALL)[3]
  cat (paste('Total Number of genes: ', nGenes, '\n', sep=''))
  cat (paste('Number of genes for downstream analysis: ', nKeep, 
             '\n', sep=''))
  expr <- expr[keepALL,,keep.lib.sizes = TRUE]
  expr = calcNormFactors(expr)
  v <- voom(expr, design=NULL, plot = FALSE)$E
  return (v)
}
