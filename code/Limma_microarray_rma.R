Limma_microarray_rma_paired <- function(eset,group, comparison,method='limma') {
  library(limma)
  group <- factor(group)
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  contrast.matrix <- makeContrasts(contrasts=comparison, 
                                   levels=design)
  # contrast.matrix<-makeContrasts(paste0(unique(group_list),
  #                                       collapse = "-"),levels = design)
  fit <- lmFit(eset, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  DEGAll <- topTable(fit2, coef=1, n = Inf)
  colnames(DEGAll) <- c('logFC', 'AveExpr', 't', 
                        'PValue', 'FDR', 'B')
  
  DEGAll$FDR <- p.adjust(DEGAll$PValue, method = 'fdr')
  o <- order(DEGAll$FDR)
  DEGAll <- DEGAll[o,]
  
  if (startsWith(rownames(eset)[1], 'ENSG')) {
    degList <- biotype[match(rownames(DEGAll), biotype$ensemblID),]
    degOutput <- data.frame(symbol=degList$geneSymbol, 
                            group=degList$group, DEGAll)
    
    keep <- which(! is.na(degOutput$symbol))
    degOutput <- degOutput[keep,]
    return(degOutput)
  } else {
    return (DEGAll)
  }
}
