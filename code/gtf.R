rm(list = ls())
options(stringsAsFactors = F)
library(rtracklayer)
gtf1 <- rtracklayer::import('gencode.v32.long_noncoding_RNAs.gtf.gz')
gtf_df1 <- as.data.frame(gtf1)
dim(gtf_df1)
save(gtf_df1,file = "gencode.v32.long_noncoding_RNAs.Rda")

gtf2 <- rtracklayer::import('gencode.v32.basic.annotation.gtf.gz')
gtf_df2 <- as.data.frame(gtf2)
dim(gtf_df2)
save(gtf_df2,file = "gencode.v32.basic.annotation.Rda")

gtf3 <- rtracklayer::import('gencode.v32.annotation.gtf.gz')
gtf_df3 <- as.data.frame(gtf3)
dim(gtf_df3)
save(gtf_df3,file = "gencode.v32.annotation.Rda")

gtf4 <- rtracklayer::import('Homo_sapiens.GRCh38.98.chr_patch_hapl_scaff.gtf.gz')
gtf_df4 <- as.data.frame(gtf4)
dim(gtf_df4)
save(gtf_df4,file = "GRCh38.98.chr_patch_hapl_scaff.Rda")

gtf5 <- rtracklayer::import('gencode.v36.annotation.gtf.gz')
gtf_df5 <- as.data.frame(gtf5)
dim(gtf_df5)
save(gtf_df5,file = "gencode.v36.annotation.Rda")
load(file="gencode.v36.annotation.Rda")


library(stringr)
genecodev36=gtf_df5
genecodev36$gene_id=str_split(genecodev36$gene_id,'[.]',simplify = T)[,1]
genecodev36$transcript_id=str_split(genecodev36$transcript_id,'[.]',simplify = T)[,1]
genecodev36$exon_id=str_split(genecodev36$exon_id,'[.]',simplify = T)[,1]
genecodev36$protein_id=str_split(genecodev36$protein_id,'[.]',simplify = T)[,1]
genecodev36$ccdsid=str_split(genecodev36$ccdsid,'[.]',simplify = T)[,1]

genecodev36_for_HRT=genecodev36[,c("gene_id","transcript_id","gene_type","gene_name")]
genecodev36_for_HRT=subset(genecodev36_for_HRT,!is.na(genecodev36_for_HRT$transcript_id))
genecodev36_for_HRT=dplyr::distinct(genecodev36_for_HRT)
save(genecodev36_for_HRT,file="genecodev36_HRT.Rdata")

genecodev36=genecodev36[,c("gene_id","gene_type","gene_name")]
genecodev36=dplyr::distinct(genecodev36)
save(genecodev36,file="genecodev36.Rdata")


library(miRBaseConverter)
library(rtracklayer)
gtf6 <- rtracklayer::import.gff3('mir_hsa_v22.gff3')
gtf_df6 <- as.data.frame(gtf6)
save(gtf_df6,file = "MIRbase.v22.annotation.Rda")
load(file="MIRbase.v22.annotation.Rda")
