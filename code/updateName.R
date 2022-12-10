library(biomaRt)
library(org.Hs.eg.db)


if(F){
all_datasets<- listDatasets(mart)
human_ensembl <- useMart(dataset="hsapiens_gene_ensembl",biomart='ENSEMBL_MART_ENSEMBL')
biomaRt_data<-getBM(attributes=c("ensembl_gene_id","hgnc_symbol","gene_biotype","external_synonym"),mart = human_ensembl)
biomaRt_data[biomaRt_data==""] <- NA
biomaRt_data=biomaRt_data[!is.na(biomaRt_data$hgnc_symbol),]
biomaRt_data=unique(biomaRt_data)
save(biomaRt_data,file="biomaRt_data.Rdata")
# listAttributes(human_ensembl)
# fuck=select(mart, keys=affy, columns=c('affy_hg_u133_plus_2','entrezgene_id'),
#        keytype='affy_hg_u133_plus_2')
#https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#how-to-build-a-biomart-query
# affyids=c("202763_at","209310_s_at","207500_at")
# getBM(attributes=c("gene_biotype", "ensembl_gene_id"),
#       filters = listFilters(mart)$name, 
#       values = affyids, 
#       mart = mart)
# listAttributes(mart)$name
}

library(AnnotationDbi)
#原来这个是alias2Symbol,属于limma包，不是特别好用，还是用biomart吧。

# checkcurrentName <- function (alias, species = "Hs", expand.symbols = FALSE) 
# {
#   alias <- as.character(alias)
#   egALIAS2EG <- tryCatch(getFromNamespace("org.Hs.egALIAS2EG", "org.Hs.eg.db"), error = function(e) FALSE)
#   egSYMBOL <- tryCatch(getFromNamespace("org.Hs.egSYMBOL", "org.Hs.eg.db"), error = function(e) FALSE)
#   
#   isSymbol <- alias %in% AnnotationDbi::Rkeys(egSYMBOL)
#   alias2 <- intersect(alias[!isSymbol], AnnotationDbi::Rkeys(egALIAS2EG))
#   eg <- AnnotationDbi::mappedLkeys(egALIAS2EG[alias2])
#   c(alias[isSymbol], AnnotationDbi::mappedRkeys(egSYMBOL[eg]))
# }

# currentName=AnnotationDbi::Rkeys(egSYMBOL)
# AliasName=AnnotationDbi::Rkeys(egALIAS2EG)
# save(currentName,AliasName,file ="checklist.Rdata")
# ALIAS=rownames(probes_expr_GSE119794_rna)
updateName <- function(ALIAS){
  # #在我们已知直接使用symbol合并效果不佳
  # table(ALIAS$ensembl_gene_id!=ALIAS$ENSEMBL,useNA="ifany")
  # #gene_newlist 其实是一个对照表。都是对的不比特别注意 org.Hs.eg.db tcga 用的这个。 biomart更新更好，但不必要刻意追求
  # table(ALIAS$SYMBOL!=ALIAS$hgnc_symbol,useNA="ifany")
  # 简而言之 使用org.Hs.eg.db更新基因名从ALIAS映射到SYMBOL，再用biomart得到最新的注释
  ALIAS=as.data.frame(ALIAS)
  load("biomaRt_data.Rdata")
  ALIAS_new=ALIAS[ALIAS$ALIAS %in% biomaRt_data$hgnc_symbol,]
  ALIAS_new=unique(ALIAS_new)
  ALIAS_new=as.data.frame(ALIAS_new)
  ALIAS_new=merge(ALIAS_new,biomaRt_data,by.x="ALIAS_new",by.y="hgnc_symbol",all.x=T)
  ALIAS_new$ALIAS=ALIAS_new$ALIAS_new
  ALIAS_new=ALIAS_new[,c("ensembl_gene_id","ALIAS_new","ALIAS","gene_biotype")]
  colnames(ALIAS_new)=c("ENSEMBL","SYMBOL","ALIAS","gene_biotype")
  #
  ALIAS_old=ALIAS[!ALIAS$ALIAS %in% biomaRt_data$hgnc_symbol,]
  ALIAS_old=as.data.frame(ALIAS_old)
  ALIAS_old=merge(ALIAS_old,biomaRt_data,by.x="ALIAS_old",by.y="external_synonym",all.x=T)
  ALIAS_old=ALIAS_old[,c("ensembl_gene_id","hgnc_symbol","ALIAS_old","gene_biotype")]
  colnames(ALIAS_old)=c("ENSEMBL","SYMBOL","ALIAS","gene_biotype")
  ALIAS=rbind(ALIAS_new,ALIAS_old)
  ALIAS=ALIAS[!is.na(ALIAS$ENSEMBL),]
  ALIAS=ALIAS[!is.na(ALIAS$SYMBOL),]
  ALIAS=ALIAS[!is.na(ALIAS$ALIAS),]
  ALIAS=ALIAS[!is.na(ALIAS$gene_biotype),]
  ALIAS=unique(ALIAS)
  return(ALIAS)
  }



# alias2SymbolUsingNCBI(alias=c("C12orf80","food","AA06"),gene.info.file = "Homo_sapiens.gene_info.gz")
# alias2SymbolUsingNCBImodify(alias=c("C12orf80","ARNT","FOOD","AA06"),gene.info.file = "Homo_sapiens.gene_info.gz")


alias2SymbolUsingNCBImodify <- function(alias,gene.info.file,required.columns=c("GeneID","Symbol","ENSEMBL","gene_type"))
{ 
  library(stringr)
  library(tidyr)
  alias <- as.character(alias)
  gene.info.file <- as.character("Homo_sapiens.gene_info.gz")
  NCBI <- read.delim(gene.info.file,comment.char="",quote="",colClasses="character")
  load("genecodev36.Rdata")
  NCBI1=subset(NCBI,grepl("Ensembl:ENSG",NCBI$dbXrefs,ignore.case=T))
  NCBI1=separate_rows(NCBI1,dbXrefs,sep = "\\|")
  NCBI1=subset(NCBI1,grepl("Ensembl:ENSG",NCBI1$dbXrefs,ignore.case=T))
  NCBI1$dbXrefs=str_split(NCBI1$dbXrefs,":",simplify = T)[,2]
  NCBI1=merge(NCBI1,genecodev36,by.x="dbXrefs",by.y="gene_id")
  NCBI1=NCBI1[,c("dbXrefs","GeneID","Symbol","Synonyms","gene_type")]
  colnames(NCBI1)=c("ENSEMBL","GeneID","Symbol","Synonyms","gene_type")
  NCBI2=subset(NCBI,!grepl("Ensembl:ENSG",NCBI$dbXrefs,ignore.case=T))
  NCBI2=merge(NCBI2,genecodev36,by.x="Symbol",by.y="gene_name")
  NCBI2=NCBI2[,c("gene_id","GeneID","Symbol","Synonyms","gene_type")]
  colnames(NCBI2)=c("ENSEMBL","GeneID","Symbol","Synonyms","gene_type")
  NCBI=rbind(NCBI1,NCBI2)
  # NCBI_genecode=separate_rows(NCBI_genecode,dbXrefs,sep = "\\|")
  # NCBI_genecode=subset(NCBI_genecode,grepl("Ensembl:ENSG",NCBI_genecode$dbXrefs,ignore.case=T)|grepl("ENSG",NCBI_genecode$gene_id,ignore.case=T))
  # NCBI_genecode1=NCBI_genecode[is.na(NCBI_genecode$gene_id),]
  # NCBI_genecode1=subset(NCBI_genecode1,grepl("Ensembl:ENSG",NCBI_genecode1$dbXrefs,ignore.case=T))
  # NCBI_genecode1$gene_id=str_split(NCBI_genecode1$dbXrefs,":",simplify = T)[,2]
  # NCBI_genecode2=NCBI_genecode[!is.na(NCBI_genecode$gene_id),]
  # NCBI_genecode=rbind(NCBI_genecode1,NCBI_genecode2)
  # table(is.na(NCBI_genecode$gene_id))
  
  #	Try matching to symbols
  m <- match(alias,NCBI$Symbol)
  EntrezID <- NCBI$GeneID[m]
  
  #	For any rows that don't match symbols, try synonyms
  i <- which(is.na(EntrezID))
  if(any(i)) {
    S <- strsplit(NCBI$Synonyms,split="\\|")
    N <- vapply(S,length,1)
    Index <- rep.int(1:nrow(NCBI),times=N)
    IS <- data.frame(Index=Index,Synonyms=unlist(S),stringsAsFactors=FALSE)
    m <- match(alias[i],IS$Synonyms)
    EntrezID[i] <- NCBI$GeneID[IS$Index[m]]
    m <- match(EntrezID,NCBI$GeneID)
  }
  
  NCBI[m,required.columns,drop=FALSE]
}
# updateName <- function(ALIAS){
#   # #在我们已知直接使用symbol合并效果不佳
#   # table(ALIAS$ensembl_gene_id!=ALIAS$ENSEMBL,useNA="ifany")
#   # #gene_newlist 其实是一个对照表。都是对的不比特别注意 org.Hs.eg.db tcga 用的这个。 biomart更新更好，但不必要刻意追求
#   # table(ALIAS$SYMBOL!=ALIAS$hgnc_symbol,useNA="ifany")
#   # 简而言之 使用org.Hs.eg.db更新基因名从ALIAS映射到SYMBOL，再用biomart得到最新的注释
#   ALIAS=as.data.frame(ALIAS)
#   # load("checklist.Rdata")
#   load("biomaRt_data.Rdata")
#   ALIAS_new=ALIAS[ALIAS$ALIAS %in% biomaRt_data$hgnc_symbol,]
#   ALIAS_new=unique(ALIAS_new)
#   # df=biomaRt::select(org.Hs.eg.db, keys = ALIAS_new, columns=c("ENSEMBL","SYMBOL"), keytype="SYMBOL")
#   ALIAS_new=as.data.frame(ALIAS_new)
#   ALIAS_new=merge(ALIAS_new,biomaRt_data,by.x="ALIAS_new",by.y="hgnc_symbol",all.x=T)
#   
#   colnames(ALIAS_new)=c("SYMBOL","ENSEMBL")
#   ALIAS_new=ALIAS_new[!is.na(ALIAS_new$SYMBOL),]
#   
#   
#   
#   # ALIAS_new[is.na(ALIAS_new$ensembl_gene_id),"ensembl_gene_id"] <- "vvv"
#   # table(is.na(biomaRt_data$ensembl_gene_id),useNA="ifany")
#   # table(is.na(biomaRt_data$hgnc_symbol),useNA="ifany")
#   # table(is.na(ALIAS_new$SYMBOL),useNA="ifany")
#   # table(is.na(ALIAS_new$ENSEMBL),useNA="ifany")#24051   691 
#   # table(is.na(ALIAS_new$ensembl_gene_id),useNA="ifany")#24051   691 
#   ALIAS_new=merge(ALIAS_new,biomaRt_data,by.x="ENSEMBL",by.y="ensembl_gene_id",all.x=T)
#   ALIAS_new$ALIAS=ALIAS_new$SYMBOL
#   ALIAS_new=ALIAS_new[,c("ENSEMBL","SYMBOL","ALIAS","gene_biotype")]
#   ALIAS_new=unique(ALIAS_new)
#   # load("biomaRt_data.Rdata")
#   # biomaRt_data=biomaRt_data[!is.na(biomaRt_data$hgnc_symbol),]
#   ALIAS_old=ALIAS[!ALIAS$ALIAS %in% currentName,]
#   ALIAS_old=unique(ALIAS_old)
#   # df=biomaRt::select(org.Hs.eg.db, keys = ALIAS_old, columns=c("ENSEMBL","SYMBOL"), keytype="ALIAS")
#   ALIAS_old=as.data.frame(ALIAS_old)
#   ALIAS_old=merge(ALIAS_old,df,by.x="ALIAS_old",by.y="ALIAS",all.x=T)
#   ALIAS_old=ALIAS_old[!is.na(ALIAS_old$SYMBOL),]
#   # ALIAS_old[is.na(ALIAS_old$ensembl_gene_id),"ensembl_gene_id"] <- "vvv"
#   # table(is.na(biomaRt_data$ensembl_gene_id),useNA="ifany")
#   # table(is.na(biomaRt_data$hgnc_symbol),useNA="ifany")
#   # table(is.na(ALIAS_old$SYMBOL),useNA="ifany")
#   # table(is.na(ALIAS_old$ENSEMBL),useNA="ifany")#24051   691 
#   # table(is.na(ALIAS_old$ensembl_gene_id),useNA="ifany")#24051   691 
#   ALIAS_old=merge(ALIAS_old,biomaRt_data,by.x="ENSEMBL",by.y="ensembl_gene_id",all.x=T)
#   #
#   ALIAS_old_1=ALIAS_old[is.na(ALIAS_old$hgnc_symbol),]
#   ALIAS_old_1=ALIAS_old_1[,-which(colnames(ALIAS_old_1) %in% c("hgnc_symbol","gene_biotype"))]
#   ALIAS_old_1=merge(ALIAS_old_1,biomaRt_data,by.x="SYMBOL",by.y="hgnc_symbol",all=F)
#   ALIAS_old_1$ENSEMBL=ALIAS_old_1$ensembl_gene_id
#   ALIAS_old_1=ALIAS_old_1[,c("ENSEMBL","SYMBOL","ALIAS_old","gene_biotype")]
#   colnames(ALIAS_old_1)=c("ENSEMBL","SYMBOL","ALIAS","gene_biotype")
#   # ensg
#   ALIAS_old_2=ALIAS_old[!is.na(ALIAS_old$hgnc_symbol),]
#   ALIAS_old_2=ALIAS_old_2[,c("ENSEMBL","SYMBOL","ALIAS_old","gene_biotype")]
#   colnames(ALIAS_old_2)=c("ENSEMBL","SYMBOL","ALIAS","gene_biotype")
#   ALIAS_old=rbind(ALIAS_old_1,ALIAS_old_2)
#   ALIAS_old=unique(ALIAS_old)
#   
#   ALIAS=rbind(ALIAS_new,ALIAS_old)
#   ALIAS=unique(ALIAS)
#   return(ALIAS)
#   
# }

if(F){
#不好用没学会
# library(biomaRt)
# library(GEOquery)
# affy=c("202763_at","209310_s_at","207500_at")
# select(mart, keys=affy, columns=c('affy_hg_u133_plus_2','entrezgene_id'),
#        keytype='affy_hg_u133_plus_2')
# 
# mart <- useMart(dataset="hsapiens_gene_ensembl",biomart='ensembl')
# 
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#how-to-build-a-biomart-query
df=biomaRt::select(org.Hs.eg.db, keys = rownames(genes_expr_mad_GSE43797_rna), columns=c("ENSEMBL","SYMBOL"), keytype="ALIAS")
# idType(OrgDb = "org.Hs.eg.db")
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
# [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"        
# [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
# [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"     
# [26] "UNIPROT" 
}

