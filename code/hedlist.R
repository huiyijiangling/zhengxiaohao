# housekeeping genes
# crispr dependency gene
# essential genes

load(file = r"(C:\Users\zxh\Desktop\x学习笔记\必需基因 管家基因\HRT Atlas v1.0\human\Housekeeping_GenesHuman.RData)")
HRT_V1_human_all=Housekeeping_Genes
Housekeeping_Genes=NULL
HRT_V1_human_most20=read.csv(file = r"(C:\Users\zxh\Desktop\x学习笔记\必需基因 管家基因\HRT Atlas v1.0\human\MostStable.csv)",sep = ";")
load(file=r"(C:\Users\zxh\Desktop\R\gtf\genecodev36_HRT.Rdata)")
table(HRT_V1_human_all$Ensembl%in%genecodev36_for_HRT$transcript_id)
HRT_V1_human_mapped=merge(HRT_V1_human_all,genecodev36_for_HRT,by.x="Ensembl",by.y="transcript_id")
table(duplicated(HRT_V1_human_mapped$Gene.name))
table(HRT_V1_human_mapped$gene_type)
HRT_V1_human_mapped=subset(HRT_V1_human_mapped,select = c("gene_id","gene_name"))
HRT_V1_human_mapped=unique(HRT_V1_human_mapped)
save(HRT_V1_human_mapped,file = "HRT_gene.Rdata")

Wang=readxl::read_xlsx(path = r"(C:\Users\zxh\Desktop\x学习笔记\必需基因 管家基因\Wang\aac7041_sm_table_s3.xlsx)")
Wang=subset(Wang,Wang$`KBM7 adjusted p-value`<0.05)
Wang=subset(Wang,Wang$`KBM7 CS`<(0))

source(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\updateName.R)")# ncbi与biomart相比多map了200，等于是，ncbi》biomart》org.hs
probe2gene=data.frame(probe_id=Wang$Gene,symbol=Wang$Gene)
forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol,gene.info.file = "Homo_sapiens.gene_info.gz")
probe2gene=cbind(probe2gene,forthetree)
Wang_ensg=subset(probe2gene,!is.na(probe2gene$ENSEMBL))
Wang_ensg=subset(Wang_ensg,select=c(Symbol,ENSEMBL))
save(Wang_ensg,file = "Wang.Rdata")


gene2ensembl_NCBI=data.table::fread(r"(C:\Users\zxh\Desktop\x学习笔记\必需基因 管家基因\单细胞用管家基因qc\gene2ensembl.gz)",header=T,na.strings = "-")
colnames(gene2ensembl_NCBI)[1]="tax_id"

# #邮箱匹配：  
# text2<-c("704232753@qq.com is my email address.")  
# grepl("[0-9.*]+@[a-z.*].[a-z.*]",text2)  
# 结果如下
# 
# > text2<-c("704232753@qq.com is my email address.")  
# > grepl("[0-9.*]+@[a-z.*].[a-z.*]",text2)  
# [1] TRUE  

gene2ensembl_NCBI=subset(gene2ensembl_NCBI,grepl("^ENSG+[0-9.*]",gene2ensembl_NCBI$Ensembl_gene_identifier,ignore.case = T))
gene2ensembl_NCBI$Ensembl_gene_identifier=str_split(gene2ensembl_NCBI$Ensembl_gene_identifier,'[.]',simplify = T)[,1]
gene2ensembl_NCBI$RNA_nucleotide_accession.version=str_split(gene2ensembl_NCBI$RNA_nucleotide_accession.version,'[.]',simplify = T)[,1]
gene2ensembl_NCBI$Ensembl_rna_identifier=str_split(gene2ensembl_NCBI$Ensembl_rna_identifier,'[.]',simplify = T)[,1]
gene2ensembl_NCBI$protein_accession.version=str_split(gene2ensembl_NCBI$protein_accession.version,'[.]',simplify = T)[,1]
gene2ensembl_NCBI$Ensembl_protein_identifier=str_split(gene2ensembl_NCBI$Ensembl_protein_identifier,'[.]',simplify = T)[,1]

gene2ensembl_NCBI_NM=subset(gene2ensembl_NCBI,grepl("^nm_",gene2ensembl_NCBI$RNA_nucleotide_accession.version,ignore.case = T))
gene2ensembl_NCBI_NM=dplyr::distinct(gene2ensembl_NCBI_NM)
gene2ensembl_NCBI_NM=subset(gene2ensembl_NCBI_NM,select=c("Ensembl_gene_identifier","RNA_nucleotide_accession.version"))
gene2ensembl_NCBI_NM=dplyr::distinct(gene2ensembl_NCBI_NM)
# gene2ensembl_NCBI$RNA_nucleotide_accession.version NM
# gene2ensembl_NCBI$Ensembl_rna_identifier ENST
Homo_HK_2013=data.table::fread(r"(C:\Users\zxh\Desktop\x学习笔记\必需基因 管家基因\单细胞用管家基因qc\2013 revisited\hg38.HouseKeepingGenes.bed.gz)",header=F)
Homo_HK_2013=subset(Homo_HK_2013,select="V4")
colnames(Homo_HK_2013)="names2015"
Homo_HK_2013=merge(Homo_HK_2013,gene2ensembl_NCBI_NM,by.x="names2015",by.y="RNA_nucleotide_accession.version",all.x=T)
# Homo_HK_2013_NA=subset(Homo_HK_2013,is.na(Homo_HK_2013$Ensembl_gene_identifier))

table(Homo_HK_2013$V4 %in% gene2ensembl_NCBI_NM$RNA_nucleotide_accession.version) #FALSE  TRUE 65  3737 
table(Homo_HK_2013_old$V2 %in% gene2ensembl_NCBI_NM$RNA_nucleotide_accession.version) #FALSE  TRUE 
table(Homo_HK_2013$Ensembl_gene_identifier %in% Homo_HK_2013_ensg$ENSEMBL) #FALSE  TRUE 65  3737 


Homo_HK_2013_old=data.table::fread(r"(C:\Users\zxh\Desktop\x学习笔记\必需基因 管家基因\单细胞用管家基因qc\2013 revisited\HK_genes.txt)",header=F)
source(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\updateName.R)")# ncbi与biomart相比多map了200，等于是，ncbi》biomart》org.hs
probe2gene=data.frame(probe_id=Homo_HK_2013_old$V1,symbol=Homo_HK_2013_old$V1)
forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol,gene.info.file = "Homo_sapiens.gene_info.gz")
probe2gene=cbind(probe2gene,forthetree)
Homo_HK_2013_old=merge(probe2gene,Homo_HK_2013_old,,by.x="probe_id",by.y="V1",all.x=T)
colnames(Homo_HK_2013_old)[ncol(Homo_HK_2013_old)]="names2013"
table(Homo_HK_2013$names2015 %in% Homo_HK_2013_old$names2013)#都一样啊
Homo_HK_2013=merge(Homo_HK_2013_old,Homo_HK_2013,by.x="names2013",by.y="names2015")
table(Homo_HK_2013$ENSEMBL==Homo_HK_2013$Ensembl_gene_identifier)#10
Homo_HK_2013_buyiyang=subset(Homo_HK_2013,(Homo_HK_2013$ENSEMBL!=Homo_HK_2013$Ensembl_gene_identifier)|is.na(Homo_HK_2013$ENSEMBL))
openxlsx::write.xlsx(Homo_HK_2013_buyiyang,file = "Homo_HK_2013手动补充.xlsx")
# 发现新版本 其实更加准备，但是其实呢，现在绝对不能更新，因为这样会导致 前后有可能叠加不上 所以要用旧版
# 注意以后要改啊。
Homo_HK_2013=subset(Homo_HK_2013,!is.na(Homo_HK_2013$ENSEMBL))
Homo_HK_2013=subset(Homo_HK_2013,select=c(Symbol,ENSEMBL))
save(Homo_HK_2013,file = "Homo_HK_2013.Rdata")
