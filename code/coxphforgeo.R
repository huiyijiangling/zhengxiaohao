coxphforgeo <- function(genes, rna.expr, metaMatrix) {
  samples = intersect(colnames(rna.expr), metaMatrix$geo_accession)
  exprDa=rna.expr[genes,samples]
  
  clinicalDa=metaMatrix[match(samples,metaMatrix$geo_accession),]
  daysToDeath <- as.numeric(clinicalDa$time)
  vitalStatus <- as.numeric(ifelse(as.numeric(clinicalDa$event)==0, 0, 1))
  exprDa=as.data.frame(exprDa)
  coxphDEGs <- c()
  for (i in seq_len(nrow(exprDa))) {
    DEG <- unlist(exprDa[i,])
    coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ DEG)
    
    summcph <- summary(coxtest)
    coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
                summcph$coefficients[,5])
    coxphDEGs <- rbind(coxphDEGs, coeffs)
    
  }
  rownames(coxphDEGs) <- rownames(exprDa)
  
  colnames(coxphDEGs) <- c('coef','HR','lower95','upper95','pValue')
  coxphDEGs <- data.frame(coxphDEGs)
  #coxphDEGs$FDR <- p.adjust(coxphDEGs$pValue, method='fdr')
  
  #o <- order(coxphDEGs$pValue)
  #coxphDEGs <- coxphDEGs[o,]
  
  return (coxphDEGs)
}
k1=dePC[mrna_2list,]
k2=genes_expr_mean_GSE84433[which( rownames(genes_expr_mean_GSE84433)%in% k1$symbol),]
survOutput <- coxphforgeo(genes=rownames(k2),
                     rna.expr=genes_expr_mean_GSE84433,
                     metaMatrix=v2)

survOutput$pValue=as.numeric(survOutput$pValue)
XXXXXXXXXXXXXXXXx=rownames(survOutput[survOutput$pValue<0.05,])

mrna_3listcox=XXXXXXXXXXXXXXXXx

exprSet_union=genes_expr_mean_GSE84433[XXXXXXXXXXXXXXXXx,]

#
# load("C:/Users/zxh/Desktop/R/gastric cancer validation/GSE62254_after_soft.Rdata")
# load("C:/Users/zxh/Desktop/R/gastric cancer validation/GSE84433_after_soft.Rdata")
# load("C:/Users/zxh/Desktop/R/gastric cancer validation/GSE29272_after_bioc.Rdata")
load("C:/Users/zxh/Desktop/R/gastric cancer validation/GSE62254_after_bioc.Rdata")
load("C:/Users/zxh/Desktop/R/gastric cancer validation/GSE84433_after_bioc.Rdata")
#check
dim(genes_expr_mean_GSE62254[which(rownames(genes_expr_mean_GSE62254) %in% rownames(exprSet_union)),])#bioc127>124
dim(genes_expr_mean_GSE84433[which(rownames(genes_expr_mean_GSE84433) %in% rownames(exprSet_union)),])#bioc127>109
# dim(genes_expr_mean_GSE29272[which(rownames(genes_expr_mean_GSE29272) %in% rownames(exprSet_union)),])#
#
symbol_union=list(exprSet_union=dePC$symbol,
                  genes_expr_mean_GSE84433=rownames(genes_expr_mean_GSE84433),
                  genes_expr_mean_GSE62254=rownames(genes_expr_mean_GSE62254)
                  # genes_expr_mean_GSE29272=rownames(genes_expr_mean_GSE29272)
)
symbol_union=Reduce(intersect,symbol_union)
exprSet_union=exprSet_union[which(rownames(exprSet_union) %in% symbol_union),]
#
genes_expr_mean_GSE62254_bioc_union=genes_expr_mean_GSE62254[which(rownames(genes_expr_mean_GSE62254) %in% rownames(exprSet_union)),]
genes_expr_mean_GSE84433_bioc_union=genes_expr_mean_GSE84433[which(rownames(genes_expr_mean_GSE84433) %in% rownames(exprSet_union)),]
# genes_expr_mean_GSE29272_bioc_union=genes_expr_mean_GSE29272[which(rownames(genes_expr_mean_GSE29272) %in% rownames(exprSet_union)),]
#
exprSet2=exprSet_union#exprSet[mrna_3listcox,]
exprSet3=as.data.frame(t(exprSet2))
exprSet4=exprSet3
colnames(exprSet4)=paste0("A",1:ncol(exprSet4))
exprSet4$geo_accession=rownames(exprSet4)
#8
dat=merge(phenoDat_GSE84433,exprSet4,by="geo_accession")

paste(colnames(exprSet4),collapse = "+")
paste(colnames(exprSet4),collapse = " ")
write.csv(dat,file="lx.csv",quote=T)
llllll=exprSet3[,c(6,47,48,51,58,66)]
llllll2=dePC[which(dePC$symbol %in% colnames(llllll)),]
paste(llllll2$symbol,collapse = "+")
paste(llllll2$symbol,collapse = "','")
