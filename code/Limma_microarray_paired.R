Limma_microarray_paired <- function(eset, group, comparison,paired,method='limma') {
  library(limma)
  group <- factor(group)
  paired <- factor(paired) 
  design <- model.matrix(~paired+group)
  # colnames(design) <- levels(group)
  # contrast.matrix <- makeContrasts(contrasts=comparison, 
  #                                  levels=design)
  # contrast.matrix<-makeContrasts(paste0(unique(group_list),
  #                                       collapse = "-"),levels = design)
  fit <- lmFit(eset, design)
  fit2 <- eBayes(fit)
  DEGAll <- topTable(fit2, coef=paste0("group",strsplit(comparison,"-")[[1]][1]), n = Inf)
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
