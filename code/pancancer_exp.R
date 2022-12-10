rm(list=ls())
options(stringsAsFactors = F)
gc()
# pancancer_exp=readr::read_tsv(file=r"(C:\Users\zxh\Desktop\R\meta collect\pcawg_publish\source\EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv)",guess_max=min(1000000, Inf))
# library_information=readr::read_tsv(file="C:/Users/zxh/Desktop/耐药文献/icsdb/library_information.tsv",guess_max=min(1000000, Inf))
# screen_information=readr::read_tsv(file="C:/Users/zxh/Desktop/耐药文献/icsdb/screen_information.tsv",guess_max=min(1000000, Inf))

#pcawg

pcawg_maf = maftools::read.maf(r"(C:\Users\zxh\Desktop\R\meta collect\pcawg_publish\source\mc3.v0.2.8.PUBLIC.maf.gz)")
pcawg_RPPA<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\meta collect\pcawg_publish\source\TCGA-RPPA-pancan-clean.txt)",header = T,verbose = F, integer64 = 'numeric')


broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg#cp
TCGA_mastercalls.abs_segtabs.fixed.txt#abs
TCGA_mastercalls.abs_tables_JSedit.fixed.txt#abs
#
system.time(pcawg_exp <- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\meta collect\pcawg_publish\source\EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv)",header = T,verbose = F, integer64 = 'numeric'))

pcawg_clinical_CDR <- readxl::read_excel(r"(C:\Users\zxh\Desktop\R\meta collect\pcawg_publish\source\TCGA-CDR-SupplementalTableS1.xlsx)")
pcawg_clinical_CDR$sampleadd=ifelse(pcawg_clinical_CDR$type=="LAML","-03","-01")
pcawg_clinical_CDR$sample=paste0(pcawg_clinical_CDR$bcr_patient_barcode,pcawg_clinical_CDR$sampleadd)
pcawg_clinical_CDR=as.data.frame(pcawg_clinical_CDR)
rownames(pcawg_clinical_CDR)=pcawg_clinical_CDR$sample
#
pcawg_clinical_PANCAN<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\meta collect\pcawg_publish\source\clinical_PANCAN_patient_with_followup.tsv)",header = T,verbose = F, integer64 = 'numeric')
#
pcawg_quality_annotations<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\meta collect\pcawg_publish\source\merged_sample_quality_annotations.tsv)",header = T,verbose = F, integer64 = 'numeric',na.strings = "")
table(pcawg_quality_annotations$Do_not_use,useNA = "ifany")
# FALSE  TRUE 
# 76672  2614 
table(pcawg_quality_annotations$AWG_excluded_because_of_pathology,useNA = "ifany")
# 0     1  <NA> 
#   62915  1179 15192 
pcawg_quality_annotations$AWG_excluded_because_of_pathology[is.na(pcawg_quality_annotations$AWG_excluded_because_of_pathology)] <- 0
table(pcawg_quality_annotations$AWG_excluded_because_of_pathology,useNA = "ifany")
pcawg_quality_annotations=subset(pcawg_quality_annotations,
                                 (is.na(pcawg_quality_annotations$sample_annotation)
                                  |(!grepl("Notification",pcawg_quality_annotations$sample_annotation)))& 
                                   (is.na(pcawg_quality_annotations$aliquot_annotation)
                                    |(!grepl("CenterNotification",pcawg_quality_annotations$aliquot_annotation))&
                                      is.na(pcawg_quality_annotations$aliquot_annotation_updated))&
                                   pcawg_quality_annotations$Do_not_use==F&
                                 (pcawg_quality_annotations$`cancer type`=="LAML"|(
                                      (is.na(pcawg_quality_annotations$patient_annotation)|(!grepl("Observation|Notification",pcawg_quality_annotations$patient_annotation)))&
                                pcawg_quality_annotations$AWG_excluded_because_of_pathology==0&
                                is.na(pcawg_quality_annotations$AWG_pathology_exclusion_reason)))
)
table(pcawg_quality_annotations$aliquot_annotation_updated)
# https://github.com/umccr/TCGA-data-harmonization/blob/master/expression/README.md
# https://github.com/umccr/TCGA-data-harmonization/blob/master/expression/TCGA_projects_summary.md
# table(pcawg_quality_annotations$Do_not_use)
aa=subset(pcawg_quality_annotations,Do_not_use==TRUE)
# cc=subset(pcawg_quality_annotations,!is.na(pcawg_quality_annotations$aliquot_annotation_updated))
dd=subset(pcawg_quality_annotations,!is.na(pcawg_quality_annotations$platform))
length(unique(dd$patient_barcode))
bb=subset(pcawg_quality_annotations,grepl("RNASeqV2",pcawg_quality_annotations$platform,ignore.case = T))
nrow(bb)#10013
nrow(subset(bb,substr(bb$aliquot_barcode,14,15)%in%c("01","03")))
# 8966
# 到 9991
# 文献9438 01，03
# table(bb$platform)
# IlluminaGA_miRNASeq    IlluminaGA_RNASeqV2 IlluminaHiSeq_miRNASeq IlluminaHiSeq_RNASeqV2 
# 1100                    557                   8118                   8820 
TCGA_tpm_ensg=list()
TCGA_tpm_trans1=list()
TCGA_tpm_trans2=list()
TCGA_tpm_symbol=list()

for(i in 1:33){
  library(dplyr)
  load(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+tpm+\\.(Rdata)",full.names = T)[i])
  # load(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+counts+\\.(Rdata)",full.names = T)[i])
  # tpm=counts
  print(ncol(tpm)-3)
  tpm=subset(tpm,select=c(T,T,T,colnames(tpm[,4:ncol(tpm)])%in%pcawg_quality_annotations$aliquot_barcode))
  print(ncol(tpm)-3)
  tpm=subset(tpm,select=c(T,T,T,!duplicated(substr(colnames(tpm[,4:ncol(tpm)]),1,15))))
  print(ncol(tpm)-3)
  tpm=subset(tpm,select=c(T,T,T,substr(colnames(tpm[,4:ncol(tpm)]),14,15)%in%c("01","03")))
  colnames(tpm)[4:ncol(tpm)]=substr(colnames(tpm[,4:ncol(tpm)]),1,15)
  tpm$means=rowMeans(tpm[,4:ncol(tpm)])
  tpm$id=stringr::str_split(tpm$id,"\\.",simplify = T)[,1]
  tpm=arrange(tpm,desc(tpm$means))
  tpm=tpm[!duplicated(tpm$id),]
  tpm=arrange(tpm,tpm$id)
  TCGA_tpm_trans1[[i]]=tpm[,1:3]
  rownames(tpm)=tpm$id
  tpm_sym=tpm[!duplicated(tpm$gene_name),]
  rownames(tpm_sym)=tpm_sym$gene_name
  tpm$means=NULL
  tpm_sym$means=NULL
  tpm=tpm[,4:ncol(tpm)]
  TCGA_tpm_trans2[[i]]=tpm_sym[,1:3]
  tpm_sym=tpm_sym[,4:ncol(tpm_sym)]
  TCGA_tpm_ensg[[i]]=tpm
  TCGA_tpm_symbol[[i]]=tpm_sym
  print(i)
}
stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+tpm+\\.(Rdata)",full.names = F),"-",simplify = T)[,2]


ensg_symbol_trans=TCGA_tpm_trans1[[1]]
save(ensg_symbol_trans,file="ensg_symbol_relation.Rdata")
names(TCGA_tpm_ensg)=stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+tpm+\\.(Rdata)",full.names = F),"-",simplify = T)[,2]
names(TCGA_tpm_symbol)=stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+tpm+\\.(Rdata)",full.names = F),"-",simplify = T)[,2]
save(pcawg_clinical_CDR,TCGA_tpm_symbol,TCGA_tpm_ensg,file="tcga_pcawg_rnaseq.Rdata")
load("tcga_pcawg_rnaseq.Rdata")

TCGA_counts_symbol=TCGA_tpm_symbol
TCGA_counts_ensg=TCGA_tpm_ensg
names(TCGA_counts_ensg)=stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+counts+\\.(Rdata)",full.names = F),"-",simplify = T)[,2]
names(TCGA_counts_symbol)=stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+counts+\\.(Rdata)",full.names = F),"-",simplify = T)[,2]
save(TCGA_counts_symbol,TCGA_counts_ensg,file="tcga_pcawg_rnaseq_counts.Rdata")

#estimate
if(F){
library(data.table)
estimate_net=list()
for (i in 1:25) {
  estimate_net[[i]]=fread(list.files(r"(C:\Users\zxh\Desktop\R\肿瘤纯度\estimate_net\)",pattern = "V2",full.names = T)[i])
}
estimate_net=rbindlist(estimate_net)
estimate_net$Purity=cos(0.6049872018+0.0001467884 * estimate_net[,3])
# TCGA_id=list()
# for(i in 1:33){
# TCGA_id[[i]]=colnames(TCGA_tpm_ensg[[i]])  
# }
TCGA_exp=do.call(cbind,TCGA_tpm_ensg)
TCGA_exp=t(TCGA_exp)
TCGA_exp=as.data.frame(TCGA_exp)
TCGA_exp$sample=stringr::str_split(rownames(TCGA_exp),"\\.",simplify = T)[,2]
TCGA_exp$project=stringr::str_split(rownames(TCGA_exp),"\\.",simplify = T)[,1]
estimate_net=merge(TCGA_exp,estimate_net,by.x="sample",by.y="ID")
rownames(estimate_net)=estimate_net$sample
save(estimate_net,file="estimate_net25.Rdata")
}
#

table(CRISPR_common_essentials$symbol%in%row.names(tpm_sym))
CRISPR_common_essentials$symbol[CRISPR_common_essentials$symbol%in%row.names(tpm_sym)]

# for(i in 1:33){
#   print(unique(pcawg_clinical_CDR[pcawg_clinical_CDR$sample%in%colnames(TCGA_tpm_ensg[[i]]),"type"])==stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+tpm+\\.(Rdata)",full.names = F),"-",simplify = T)[,2][i]
# )
#   print(unique(pcawg_clinical_CDR[pcawg_clinical_CDR$sample%in%colnames(TCGA_tpm_symbol[[i]]),"type"])==stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+tpm+\\.(Rdata)",full.names = F),"-",simplify = T)[,2][i]
# )
#   print(i)
# }
# names(TCGA_tpm_ensg)=stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+tpm+\\.(Rdata)",full.names = F),"-",simplify = T)[,2]
# names(TCGA_tpm_symbol)=stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^TCGA\\-+[[:alpha:]]+\\-+tpm+\\.(Rdata)",full.names = F),"-",simplify = T)[,2]

#cptac 2/3
#还不会写
for(i in 1:2){
  library(dplyr)
  load(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "^CPTAC\\-+[0-9]+\\-+tpm+\\.(Rdata)",full.names = T)[i])
  print(ncol(tpm)-3)
  tpm=subset(tpm,select=c(T,T,T,colnames(tpm[,4:ncol(tpm)])%in%pcawg_quality_annotations$aliquot_barcode))
  print(ncol(tpm)-3)
  tpm=subset(tpm,select=c(T,T,T,!duplicated(substr(colnames(tpm[,4:ncol(tpm)]),1,15))))
  print(ncol(tpm)-3)
  tpm=subset(tpm,select=c(T,T,T,substr(colnames(tpm[,4:ncol(tpm)]),14,15)%in%c("01","03")))
  colnames(tpm)[4:ncol(tpm)]=substr(colnames(tpm[,4:ncol(tpm)]),1,15)
  tpm$means=rowMeans(tpm[,4:ncol(tpm)])
  tpm$id=stringr::str_split(tpm$id,"\\.",simplify = T)[,1]
  tpm=arrange(tpm,desc(tpm$means))
  tpm=tpm[!duplicated(tpm$id),]
  tpm=arrange(tpm,tpm$id)
  CPTAC_tpm_trans1[[i]]=tpm[,1:3]
  rownames(tpm)=tpm$id
  tpm_sym=tpm[!duplicated(tpm$gene_name),]
  rownames(tpm_sym)=tpm_sym$gene_name
  tpm$means=NULL
  tpm_sym$means=NULL
  tpm=tpm[,4:ncol(tpm)]
  CPTAC_tpm_trans2[[i]]=tpm_sym[,1:3]
  tpm_sym=tpm_sym[,4:ncol(tpm)]
  CPTAC_tpm_ensg[[i]]=tpm
  CPTAC_tpm_symbol[[i]]=tpm_sym
  print(i)
}






# FALSE  TRUE 
# 5  1907 
#先用count
#再用tpm
#应该没有差异
#
#直接引用处理成生存意义
#下载tpm
#crispr
#应该不是common essentials 这个是和rnai共有的

list.files(r"(C:\Users\zxh\Desktop\R\depmap\22Q2\)")
# CRISPR_gene_firstpage<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\DOWNLOAD FROM FIRST PAGE\CRISPR_(DepMap_22Q2_Public+Score,_Chronos).csv)",header = T,verbose = F, integer64 = 'numeric')
# first page = meta + gene effect
# dependency 超过0.5是更加可能的
# common_essentials_pos_control=data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\common_essentials.csv)",header = T,verbose = F, integer64 = 'numeric')
CRISPR_gene_dependency<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\CRISPR_gene_dependency.csv)",header = T,verbose = F, integer64 = 'numeric')
CRISPR_gene_effect<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\CRISPR_gene_effect.csv)",header = T,verbose = F, integer64 = 'numeric')
CRISPR_gene_effect=as.data.frame(CRISPR_gene_effect)
rownames(CRISPR_gene_effect)=CRISPR_gene_effect$DepMap_ID
# colnames(CRISPR_gene_effect)%in%
coe=subset(CRISPR_gene_effect,select=colnames(CRISPR_gene_effect)%in%CRISPR_common_essentials$gene)
coe
psych::describe(coe)
class(coe$`AAMP (14)`)
heatmap(as.matrix(coe))
surv_res[[1]]

load(r"(C:\Users\zxh\Desktop\R\paad-tcga-gtex\pcawg_surv_res.Rdata)")
pcawg_surv_res=surv_res
surv_res_map <- list()
for (x in 1:33) {
surv_res_map[[x]] <- as.data.frame(t(pcawg_surv_res[[x]]))
surv_res_map[[x]]=tibble::remove_rownames(surv_res_map[[x]])
}
surv_res_map=data.table::rbindlist(surv_res_map)
surv_res_map=as.data.frame(t(surv_res_map))

HR_map=surv_res_map[,seq(from=4,length.out=33,by=7)]
HR_map=as.data.frame(lapply(HR_map,as.numeric))
colnames(HR_map)=stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "tpm",full.names = F),"-",simplify = T)[,2]
rownames(HR_map)=surv_res_map$V1
eec1=na.omit(HR_map)

pvalue_map=surv_res_map[,seq(from=7,length.out=33,by=7)]
pvalue_map=as.data.frame(lapply(pvalue_map,as.numeric))
colnames(pvalue_map)=stringr::str_split(list.files(r"(C:\Users\zxh\Desktop\R\meta collect\result\)",pattern = "tpm",full.names = F),"-",simplify = T)[,2]
rownames(pvalue_map)=surv_res_map$V1
eec2=na.omit(pvalue_map)



eec1[eec2>=0.05] <- NA
eec1[eec1>=10] <- 10
eec1[eec1<=0.1] <- 0.1
eec1=log10(eec1)
eec1_RowsAllNa_removed =  eec1[apply(eec1, 1, function(y) any(!is.na(y))),]
quantile(eec1_RowsAllNa_removed$ACC,0.99,na.rm = T)

uu=eec1_RowsAllNa_removed[rownames(eec1_RowsAllNa_removed)%in%coe$id,]
coe=ensg_symbol_trans[ensg_symbol_trans$gene_name%in%CRISPR_common_essentials$symbol,]
pdf(file = "aaa.pdf",height = 4,width = 4)
ComplexHeatmap::Heatmap(as.matrix(eec1_RowsAllNa_removed[rownames(eec1_RowsAllNa_removed)%in%coe$id,]),na_col = "white",cluster_rows = F,show_row_names = F)
                        dev.off()
CRISPR_common_essentials<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\CRISPR_common_essentials.csv)",header = T,verbose = F, integer64 = 'numeric',sep = ",")
CRISPR_common_essentials$symbol=stringr::str_split(CRISPR_common_essentials$gene," ",simplify = T)[,1]
CRISPR_common_essentials$geneid=stringr::str_split(CRISPR_common_essentials$gene," ",simplify = T)[,2]




# ### CCLE_expression_full.csv
# 
# Pipeline: Expression
# 
# RNAseq TPM gene expression data for all genes using RSEM. Log2 transformed, using a pseudo-count of 1.                                                                                                                                                                                                                                                                           
 
# ### CCLE_RNAseq_reads.csv
# 
# Pipeline: Expression
# 
# RNAseq read count data from RSEM.


# CCLE_wes_gene_cn 基因少细胞系多
CCLE_gene_cn <- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\CCLE_gene_cn.csv)",header = T,verbose = F, integer64 = 'numeric')
# CCLE_wes_gene_cn 基因多细胞系少
CCLE_wes_gene_cn<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\CCLE_wes_gene_cn.csv)",header = T,verbose = F, integer64 = 'numeric')
# CCLE_wes_segment_cn 不全
CCLE_segment_cn<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\CCLE_segment_cn.csv)",header = T,verbose = F, integer64 = 'numeric')
# CCLE_wes_segment_cn 不全
# CCLE_wes_segment_cn<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\CCLE_wes_segment_cn.csv)",header = T,verbose = F, integer64 = 'numeric')
sample_info<- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\sample_info.csv)",header = T,verbose = F, integer64 = 'numeric')



# load("iCSDB_PANC_top25.Rdata")

library(tidyr)
#NC
if(F){
  SSDs=readxl::read_excel("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/StronglySelectiveDependencies ENSG.xlsx")
  common983=readxl::read_excel("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/983.xlsx")
  wuyusile=common983[common983$ENSEMBL %in% mhighall,]
}
######################################################################
if(T){
  library(stringr)
  cell_line_information=readr::read_tsv(file="C:/Users/zxh/Desktop/耐药文献/icsdb/cell_line_information.tsv",guess_max=min(1000000, Inf))
  library_information=readr::read_tsv(file="C:/Users/zxh/Desktop/耐药文献/icsdb/library_information.tsv",guess_max=min(1000000, Inf))
  screen_information=readr::read_tsv(file="C:/Users/zxh/Desktop/耐药文献/icsdb/screen_information.tsv",guess_max=min(1000000, Inf))
  screen_information=as.data.frame(screen_information)
  library_information=as.data.frame(library_information)
  cell_line_information=as.data.frame(cell_line_information)
  library_information[2,1] <- "lib1"
  filenames_1 <- list.files('C:/Users/zxh/Desktop/耐药文献/icsdb/ewid.data_gene.PubMed_bcscore/',pattern = ".csv",full.name=TRUE)
  filenames_1
  resultALL1 <- lapply(filenames_1, function(fl) 
    read.delim(fl,header = T,row.names = 1))
  for (i in 1:length(resultALL1)){
    names(resultALL1)[i] <- str_split(str_split(filenames_1,"/",simplify = T)[,8],"\\.",simplify = T)[,1][i]}
  ewid.data_gene.PubMed=resultALL1
  #
  filenames_1 <- list.files('C:/Users/zxh/Desktop/耐药文献/icsdb/ewid.data_gene.GeCKO_bcscore/',pattern = ".csv",full.name=TRUE)
  filenames_1
  resultALL1 <- lapply(filenames_1, function(fl) 
    read.delim(fl,header = T,row.names = 1))
  for (i in 1:length(resultALL1)){
    names(resultALL1)[i] <- str_split(str_split(filenames_1,"/",simplify = T)[,8],"\\.",simplify = T)[,1][i]}
  ewid.data_gene.GeCKO_bcscore=resultALL1
  #
  filenames_1 <- list.files('C:/Users/zxh/Desktop/耐药文献/icsdb/ewid.data_gene.DepMap_20Q2_bcscore/',pattern = ".csv",full.name=TRUE)
  filenames_1
  resultALL1 <- lapply(filenames_1, function(fl) 
    read.delim(fl,header = T,row.names = 1))
  for (i in 1:length(resultALL1)){
    names(resultALL1)[i] <- str_split(str_split(filenames_1,"/",simplify = T)[,8],"\\.",simplify = T)[,1][i]}
  ewid.data_gene.DepMap_20Q2_bcscore=resultALL1
  #
  filenames_1 <- list.files('C:/Users/zxh/Desktop/耐药文献/icsdb/ewid.data_gene.Sanger_bcscore/',pattern = ".csv",full.name=TRUE)
  filenames_1
  resultALL1 <- lapply(filenames_1, function(fl) 
    read.delim(fl,header = T,row.names = 1))
  for (i in 1:length(resultALL1)){
    names(resultALL1)[i] <- str_split(str_split(filenames_1,"/",simplify = T)[,8],"\\.",simplify = T)[,1][i]}
  ewid.data_gene.Sanger_bcscore=resultALL1
}
###################################################################
if(F){
  # ASPC1、SUIT2、SU8686、LS513
  # cell_line_information
  # library_information
  # screen_information
  # 
  # ewid.data_gene.Sanger_bcscore
  # ewid.data_gene.DepMap_20Q2_bcscore
  # ewid.data_gene.GeCKO_bcscore
  # ewid.data_gene.PubMed
  cell_line_information_pancreas=subset(cell_line_information,cell_line_information$CELL_LINE%in%c("ASPC1","SUIT2","SU8686","LS513"))
  screen_information_pancreas=subset(screen_information,screen_information$CELL_LINE %in%cell_line_information_pancreas$CELL_LINE)
  
  ewid.data_gene.Sanger_bcscore_1=ewid.data_gene.Sanger_bcscore [names(ewid.data_gene.Sanger_bcscore)%in%screen_information_pancreas$EWHA_ID]
  ewid.data_gene.DepMap_20Q2_bcscore_1=ewid.data_gene.DepMap_20Q2_bcscore [names(ewid.data_gene.DepMap_20Q2_bcscore)%in%screen_information_pancreas$EWHA_ID]
  ewid.data_gene.GeCKO_bcscore_1=ewid.data_gene.GeCKO_bcscore [names(ewid.data_gene.GeCKO_bcscore)%in%screen_information_pancreas$EWHA_ID]
  ewid.data_gene.PubMed_1=ewid.data_gene.PubMed [names(ewid.data_gene.PubMed)%in%screen_information_pancreas$EWHA_ID]
  
  common_list=list(ewid.data_gene.Sanger_bcscore_1,ewid.data_gene.DepMap_20Q2_bcscore_1,ewid.data_gene.GeCKO_bcscore_1,ewid.data_gene.PubMed_1)
  common_list_1=lapply(common_list,function(x) lapply(x, function(x) subset(x,!is.na(x$RANK))))
  resultTlist=unlist(common_list_1,recursive = F)
  resultTlist=lapply(resultTlist,function(x) rownames(x)[1:quantile(x$RANK,0.10000)])
  
  library(VennDiagram)#只能有五个
  # resultTlist_U1=calculate.overlap(resultTlist)
  resultTlist_U=get.venn.partitions(resultTlist)
  resultTlist_U[,c(1:length(resultTlist))]=apply(resultTlist_U[,c(1:length(resultTlist))],2,function(x) as.numeric(x))
  resultTlist_U$sum=apply(resultTlist_U[,c(1:15)],1,sum)
  ansT2fenzhi3 <- Reduce(union,resultTlist_U[resultTlist_U$sum>=(length(resultTlist)*2/3),]$..values..)
  #iCSDB  iCSDB: an integrated database of CRISPR screens  其实这里是大于等于
  # ansT4_equal <- Reduce(union,resultTlist_U[resultTlist_U$sum>3,]$..values..)
  # resultTlist_U=resultTlist_U[resultTlist_U$TCGA=="1",]
  # resultTlist_U=resultTlist_U[resultTlist_U$sum>2,]
  # ansT3 <- Reduce(union,resultTlist_U$..values..)
  # ansT4 <- Reduce(union,resultTlist_U[resultTlist_U$sum>3,]$..values..)
  # mT=stringr::str_split(ansT3,"_",simplify = T)[,2]
  # lncT=stringr::str_split(ansT3,"_",simplify = T)[,1]
  # resultT_ce=subset(resultT1[[1]],resultT1[[1]]$Genes%in%mT &resultT1[[1]]$lncRNAs %in%lncT)
  # edges <- gdcExportNetwork(ceNetwork = resultT_ce, net = 'edges')
  # nodes <- gdcExportNetwork(ceNetwork = resultT_ce, net = 'nodes')
  # write.csv(edges, file='edgesT.csv',quote = F) ### Network of Cytoscape
  # write.csv(nodes, file='nodesT.csv',quote = F) ### Table of Cytoscape
}
#
cell_line_information_pancreas=subset(cell_line_information,cell_line_information$TISSUE=="PANCREAS")
# cell_line_information_pancreas$nameupper=str_split(cell_line_information_pancreas$CCLE_NAME,"_",simplify = T)[,1]
screen_information_pancreas=subset(screen_information,screen_information$CELL_LINE %in%cell_line_information_pancreas$CELL_LINE)
ewid.data_gene.Sanger_bcscore_1=ewid.data_gene.Sanger_bcscore [names(ewid.data_gene.Sanger_bcscore)%in%screen_information_pancreas$EWHA_ID]
ewid.data_gene.DepMap_20Q2_bcscore_1=ewid.data_gene.DepMap_20Q2_bcscore [names(ewid.data_gene.DepMap_20Q2_bcscore)%in%screen_information_pancreas$EWHA_ID]
ewid.data_gene.GeCKO_bcscore_1=ewid.data_gene.GeCKO_bcscore [names(ewid.data_gene.GeCKO_bcscore)%in%screen_information_pancreas$EWHA_ID]
ewid.data_gene.PubMed_1=ewid.data_gene.PubMed [names(ewid.data_gene.PubMed)%in%screen_information_pancreas$EWHA_ID]

common_list=list(ewid.data_gene.Sanger_bcscore_1,ewid.data_gene.DepMap_20Q2_bcscore_1,ewid.data_gene.GeCKO_bcscore_1,ewid.data_gene.PubMed_1)
common_list_1=lapply(common_list,function(x) lapply(x, function(x) subset(x,!is.na(x$RANK))))
resultTlist=unlist(common_list_1,recursive = F)
resultTlist=resultTlist[!names(resultTlist)%in%"ew140"]
resultTlist=resultTlist[lapply(resultTlist,function(x) dim(x)[1])>10000]
common_intersect=lapply(resultTlist,function(x) rownames(x[x[,3]<0,]))#这里选择百分
common_intersect=Reduce(intersect,common_intersect)
resultTlist=lapply(resultTlist,function(x) rownames(x)[1:quantile(x$RANK,0.10000)])#这里选择百分比,我这里选择了25%，因为10%只有62个基因
resultTlist=lapply(resultTlist, function(x) x[x%in%common_intersect])
library(VennDiagram)#只能有五个
# resultTlist_U1=calculate.overlap(resultTlist)
library(foreach)
library(doParallel)
library(readxl)
library(bigmemory)
#按照 KRAS central
if(T){
  library(data.table)
  data.table::setDTthreads(threads =4 )#, percent = 99
  
  CCLE_mutations_21Q4 <- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\21Q4\CCLE_mutations.csv)",header = T,sep = ",",verbose = F, integer64 = 'numeric')
  CCLE_mutations_21Q4$Protein_Change_simplify=stringr::str_split(CCLE_mutations_21Q4$Protein_Change,"\\.",simplify = T)[,2]
  CCLE_mutations_21Q4$cDNA_Change_simplify=stringr::str_split(CCLE_mutations_21Q4$Protein_Change,"\\.",simplify = T,)[,2]
  CCLE_mutations_21Q4[CCLE_mutations_21Q4==""] <- NA
  table(CCLE_mutations_21Q4$Variant_annotation)#damaging other conserving other non-conserving silent 
  
  gene_cn_21Q4 <- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\21Q4\CCLE_gene_cn.csv)",header = T,sep = ",",verbose = F, integer64 = 'numeric')
  gene_cn_21Q4 = as.data.frame(gene_cn_21Q4)
  rownames(gene_cn_21Q4) = gene_cn_21Q4$V1
  
  sample_info_21Q4 <- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\21Q4\sample_info.csv)",header = T,sep = ",",verbose = F, integer64 = 'numeric')
  sample_info_21Q4 = as.data.frame(sample_info_21Q4)
  sample_info_21Q4[sample_info_21Q4=="Unknown"] <- NA
  sample_info_21Q4[sample_info_21Q4=="unknown"] <- NA
  sample_info_21Q4[sample_info_21Q4==""] <- NA
  rownames(sample_info_21Q4)=sample_info_21Q4$DepMap_ID
  table(is.na(sample_info_21Q4$sample_collection_site))
  table(is.na(sample_info_21Q4$lineage))#44
  table(is.na(sample_info_21Q4$primary_disease))#44
  table(is.na(sample_info_21Q4$Subtype))#44
  table(is.na(sample_info_21Q4$lineage_subtype))#44
  ######################################################pancreas
  sample_info_21Q4_pancreas=subset(sample_info_21Q4,lineage %in% "pancreas")
  #
  gene_cn_21Q4_pancreas=subset(gene_cn_21Q4,V1%in%sample_info_21Q4_pancreas$DepMap_ID)
  gene_cn_21Q4_pancreas_driver41=gene_cn_21Q4_pancreas[,c("KRAS (3845)","TP53 (7157)","SMAD4 (4089)","CDKN2A (1029)","GATA6 (2627)")]
  gene_cn_21Q4_pancreas_driver41_dda=lapply(gene_cn_21Q4_pancreas_driver41,function(x) cut(x,breaks = c(-Inf,0.65,1.4,Inf),right =T,ordered_result = T,labels = c("Deletion","Deploid","Amplification")))
  gene_cn_21Q4_pancreas_driver41_dda=as.data.frame(gene_cn_21Q4_pancreas_driver41_dda)
  rownames(gene_cn_21Q4_pancreas_driver41_dda)=rownames(gene_cn_21Q4_pancreas_driver41)
  rownames(sample_info_21Q4_pancreas)==rownames(gene_cn_21Q4_pancreas_driver41_dda)
  #
  CCLE_mutations_21Q4_pancreas=subset(CCLE_mutations_21Q4,DepMap_ID%in%sample_info_21Q4_pancreas$DepMap_ID)
  CCLE_mutations_21Q4_pancreas_kras=subset(CCLE_mutations_21Q4_pancreas,Hugo_Symbol=="KRAS")
  CCLE_mutations_21Q4_pancreas_kras=subset(CCLE_mutations_21Q4_pancreas_kras,!Variant_annotation=="silent")
  CCLE_mutations_21Q4_pancreas_kras_protein=CCLE_mutations_21Q4_pancreas_kras[,c("DepMap_ID","Protein_Change_simplify")]
  #
  CCLE_mutations_21Q4_pancreas_kras_protein=aggregate(CCLE_mutations_21Q4_pancreas_kras_protein$Protein_Change, list(CCLE_mutations_21Q4_pancreas_kras_protein$DepMap_ID), paste0, collapse=",")
  colnames(CCLE_mutations_21Q4_pancreas_kras_protein)=c("DepMap_ID","Protein_Change_simplify")
  # gene_cn_21Q4$V1 = NULL
  sample_info_21Q4_pancreas=merge(sample_info_21Q4_pancreas,gene_cn_21Q4_pancreas_driver41_dda,by=0)#"row.names")
  sample_info_21Q4_pancreas=sample_info_21Q4_pancreas[-1]
  
  #check 有时候depmap没有及时更新导致两个样本的mutation信息没有及时删除
  sample_info_21Q4_pancreas$DepMap_ID[!sample_info_21Q4_pancreas$DepMap_ID %in%CCLE_mutations_21Q4$DepMap_ID]#[1] "ACH-001101" "ACH-001171"
  cell_line_information_pancreas$CCLE_NAME[!cell_line_information_pancreas$CCLE_NAME %in%sample_info_21Q4_pancreas$CCLE_Name]#character(0)
  #https://web.expasy.org/cellosaurus/CVCL_3567
  #https://web.expasy.org/cellosaurus/CVCL_5146
  sample_info_21Q4_pancreas=merge(sample_info_21Q4_pancreas,all.x=T,CCLE_mutations_21Q4_pancreas_kras_protein,by = "DepMap_ID")
  sample_info_21Q4_pancreas=subset(sample_info_21Q4_pancreas,!DepMap_ID %in%  c("ACH-001101","ACH-001171"))
  # sample_info_21Q4_pancreas=subset(sample_info_21Q4_pancreas,DepMap_ID %in% sample_info_21Q4$DepMap_ID)#一般性版本
  table(is.na(cell_line_information$CELL_LINE))#CELL_LINE
  sample_info_21Q4[is.na(sample_info_21Q4$stripped_cell_line_name),]#我觉得没意义了，就用stripped_cell_line_name
  table(is.na(sample_info_21Q4$stripped_cell_line_name))#CELL_LINE
  #
  cell_line_information_pancreas=merge(cell_line_information_pancreas,sample_info_21Q4_pancreas,by.x = "CELL_LINE",by.y = "stripped_cell_line_name")
  #
  duplicated(colnames(cell_line_information))
  huluwa=sample_info_21Q4_pancreas[is.na(sample_info_21Q4_pancreas$Protein_Change_simplify),]
  
  # 开始处理了按照 kras amp +/- kras 类型
  table(cell_line_information_pancreas$Protein_Change_simplify,cell_line_information_pancreas$KRAS..3845.,useNA = "ifany")
}
#正常计算不影响下方
if(T){
  library(maftools)
  laml.plus.gistic2=sample_info_21Q4[,c("DepMap_ID","stripped_cell_line_name")]
  laml.plus.gistic2=merge(CCLE_mutations_21Q4,laml.plus.gistic2,by="DepMap_ID")
  laml.plus.gistic2=subset(laml.plus.gistic2,laml.plus.gistic2$DepMap_ID %in% cell_line_information_pancreas$DepMap_ID)
  colnames(laml.plus.gistic2)[which(colnames(laml.plus.gistic2) == "Alternate_Allele")] = "Tumor_Seq_Allele2"
  colnames(laml.plus.gistic2)[which(colnames(laml.plus.gistic2) == "stripped_cell_line_name")] = "Tumor_Sample_Barcode"
  
  
  laml.clin=cell_line_information_pancreas
  laml.clin[is.na(laml.clin$Protein_Change_simplify),"Protein_Change_simplify"] <- "No"
  colnames(laml.clin)[which(colnames(laml.clin) == "CELL_LINE")] = "Tumor_Sample_Barcode"
  colnames(laml.clin)[which(colnames(laml.clin) == "KRAS..3845.")] = "KRAS_copy_number"
  colnames(laml.clin)[which(colnames(laml.clin) == "Protein_Change_simplify")] = "KRAS_Protein_Change"
  # laml.clin$`KRAS copy number`
  # paste(colnames(laml.clin),collapse = '","')
  laml.plus.gistic2 = maftools::read.maf(laml.plus.gistic2,clinicalData = laml.clin)
  oncoplot(laml.plus.gistic2, top = 5)
  
  getClinicalData(x = laml.plus.gistic2)
  oncoplot(maf = laml.plus.gistic2, top = 5)
  pdf("maf_ccle_total.pdf")
  oncoplot(
    maf = laml.plus.gistic2,
    top = 20,
    minMut = NULL,
    genes = NULL,
    altered = FALSE,
    drawRowBar = TRUE,
    drawColBar = TRUE,
    # leftBarData = aml_genes_vaf,
    # leftBarLims = c(0, 100),
    # rightBarData = laml.mutsig[,.(gene, q)],
    rightBarLims = NULL,
    topBarData = NULL,
    logColBar = FALSE,
    includeColBarCN = TRUE,
    clinicalFeatures = c("KRAS_copy_number","KRAS_Protein_Change"),#"sex","source","Achilles_n_replicates","cell_line_NNMD","culture_type","cas9_activity","sample_collection_site","primary_or_metastasis","age","lineage_subtype"),
    annotationColor = NULL,
    annotationDat = NULL,
    # pathways = pathways,
    path_order = NULL,
    selectedPathways = NULL,
    pwLineCol = "#535c68",
    pwLineWd = 1,
    draw_titv = TRUE,
    titv_col = NULL,
    showTumorSampleBarcodes = T,
    barcode_mar = 6,#space
    barcodeSrt = 90,
    gene_mar = 5,
    anno_height = 1,
    legend_height = 4,
    sortByAnnotation = TRUE,
    groupAnnotationBySize = TRUE,
    annotationOrder = NULL,
    sortByMutation = FALSE,
    keepGeneOrder = FALSE,
    GeneOrderSort = TRUE,
    sampleOrder = NULL,
    # additionalFeature = c("Tumor_Seq_Allele2", "C"),#就是标记特征的sample
    additionalFeaturePch = 20,
    additionalFeatureCol = "gray70",
    additionalFeatureCex = 0.9,
    genesToIgnore = NULL,
    removeNonMutated = TRUE,
    fill = TRUE,
    cohortSize = NULL,
    colors = NULL,
    cBioPortal = FALSE,
    bgCol = "#CCCCCC",
    borderCol = "white",
    annoBorderCol = NA,
    numericAnnoCol = NULL,#'YlOrBr' ,
    drawBox = FALSE,
    fontSize = 0.8,
    SampleNamefontSize = 1,
    titleFontSize = 1.5,
    legendFontSize = 1.2,
    annotationFontSize = 1.2,
    sepwd_genes = 0.5,
    sepwd_samples = 0.25,
    writeMatrix = FALSE,
    colbar_pathway = FALSE,
    showTitle = TRUE,
    titleText = NULL,
    showPct = TRUE
  )
  dev.off()
  pdf("maf_ccle_summary.pdf")
  plotmafSummary(laml.plus.gistic2)
  dev.off()
}
#如果直接跳到这里直接按照
if(length(resultTlist)<=32){
  resultTlist_U=get.venn.partitions(resultTlist,keep.elements = F)
  resultTlist_U[,c(1:length(resultTlist))]=apply(resultTlist_U[,c(1:length(resultTlist))],2,function(x) as.numeric(x))
  resultTlist_U$sum=apply(resultTlist_U[,c(1:15)],1,sum)
  ansT2fenzhi3 <- Reduce(union,resultTlist_U[resultTlist_U$sum>=(length(resultTlist)*2/3),]$..values..)
}
if(length(resultTlist)>32){
  resultTlist_U=unique(unlist(resultTlist))
  cl <- makeCluster(30)
  registerDoParallel(cl)
  resultTlist_matrix <- foreach(i=resultTlist_U,.combine = rbind)%dopar% lapply(i,function(x)grepl(x,resultTlist))
  stopCluster(cl)#.combine=rbind
  resultTlist_matrix=as.data.frame(resultTlist_matrix)
  rownames(resultTlist_matrix)=resultTlist_U
  # table(resultTlist_matrix[2,])
  save(resultTlist_matrix,file="resultTlist_matrix.Rdata")
  table(screen_information_pancreas$CELL_LINE %in% cell_line_information_pancreas$CELL_LINE)
  select_criteria=screen_information_pancreas[match(names(resultTlist),screen_information_pancreas$EWHA_ID),]
  select_criteria=merge(select_criteria,cell_line_information_pancreas,by.x ="CELL_LINE",by.y="CELL_LINE")
  rownames(select_criteria)=select_criteria$EWHA_ID
  select_criteria=select_criteria[names(resultTlist),]
  select_criteria=unite(select_criteria,col="subselect","Protein_Change_simplify","KRAS..3845.",sep = "_",remove = F)
  select_criteria$subselect_num=as.numeric(as.factor(select_criteria$subselect))
  table(select_criteria$subselect_num)
  table(select_criteria$subselect,useNA = "ifany")
  #
  kras_mut_cna_matrix=rbindlist(list(resultTlist_matrix$V1))
  kras_mut_cna_matrix = lapply(kras_mut_cna_matrix, function(x) {
    test1 = table(select_criteria$subselect_num[x], useNA = "always")
    data.frame(matrix(
      as.numeric(test1),
      nrow = 1,
      dimnames = list(NULL, paste0("S", names(test1)))
    ))
  })
  kras_mut_cna_matrix=data.table::rbindlist(kras_mut_cna_matrix,use.names = T,fill = T)
  kras_mut_cna_matrix=as.data.frame(kras_mut_cna_matrix)
  rownames(kras_mut_cna_matrix)=rownames(resultTlist_matrix)
  library(dplyr)
  Corrresponding_count = select_criteria %>%
    group_by(subselect_num, subselect) %>%
    summarise(Counts = n())
  Corrresponding_count
  Corrresponding_count$subselect_num_S=paste0("S",Corrresponding_count$subselect_num)
  kras_mut_cna_matrix=subset(kras_mut_cna_matrix,select=Corrresponding_count$subselect_num_S)
  colnames(kras_mut_cna_matrix)%in%select_criteria$subselect_num
  unique(select_criteria[,c("subselect","subselect_num")])
  kras_mut_cna_matrix_check=as.data.frame(apply(kras_mut_cna_matrix,1,function(x) table(x>=(Corrresponding_count$Counts*2/3))["TRUE"]==10))
  kras_mut_cna_matrix_check_T=na.omit(kras_mut_cna_matrix_check)
  kras_mut_cna_matrix_check_T=subset(kras_mut_cna_matrix_check_T,kras_mut_cna_matrix_check_T$`apply(kras_mut_cna_matrix, 1, function(x) table(x >= (Corrresponding_count$Counts * 2/3))["TRUE"] == 10)`==TRUE)
  ansT2fenzhi3=rownames(kras_mut_cna_matrix_check_T)
  #                 subselect subselect_num
  # ew1018         NA_Deploid             9
  # ew1037 G12V_Amplification             7
  # ew1038       G12C_Deploid             2
  # ew1045       G12D_Deploid             4
  # ew1143 G12D_Amplification             3
  # ew920        G12V_Deploid             8
  # ew936        G12R_Deploid             6
  # ew979  G12R_Amplification             5
  # ew188        Q61H_Deploid            10
  # ew777        G12A_Deploid             1
  # KRAS status Deletion Deploid Amplification  Type
  # G12A        0       1             0
  # G12C        0       1             0
  # G12D        0      11             7
  # G12R        0       2             2
  # G12V        0       8             6
  # Q61H        0       2             0
  # <NA>        0       3             0
  
  
  #直接筛选
  if(F){
    ansT2fenzhi3=rownames(resultTlist_matrix)[apply(resultTlist_matrix,1,function(x) table(x)["TRUE"]>50)]
    ansT2fenzhi3=rownames(resultTlist_matrix)[apply(resultTlist_matrix,1,function(x) table(x)["TRUE"]>=69)]
  }
  save(ansT2fenzhi3,file="iCSDB_PANC_top10.Rdata")
  # load(file="iCSDB_PANC_top10.Rdata")
  # save(ansT2fenzhi3,file="iCSDB_PANC_top25.Rdata")
  # save(ansT2fenzhi3,file="iCSDB_PANC_top5.Rdata")
  # parLapply(cl = NULL, X, fun, ..., chunk.size = NULL)
  # system.time({a = 1+1})
}
load("iCSDB_PANC_top10.Rdata")
dependency_heat=unlist(common_list_1,recursive = F)
dependency_heat=dependency_heat[!names(dependency_heat)%in%"ew140"]
dependency_heat=lapply(dependency_heat,function(x) x[3])
dependency_heat=dependency_heat[lapply(dependency_heat, function(x) lengths(x))>10000]
dependency_heat=lapply(dependency_heat,function(x) x[ansT2fenzhi3,])
dependency_heat=as.data.frame(dependency_heat)
table(dependency_heat$ew1038<0)
table(is.na(dependency_heat$ew35))
rownames(dependency_heat)=ansT2fenzhi3
library(pheatmap)
library(viridis)
pdf("oridinal no ew140 plot1.pdf")
pheatmap(
  mat               = dependency_heat,#t(cell_tumor_matrix)
  color             = viridis(50),#inferno(50),
  # border_color      = NA,
  show_colnames     = T,
  show_rownames     = F,
  fontsize          = 8,
  cluster_cols = T,
  # annotation_col = annotation_col,
  annotation_legend = T,
  # annotation_colors = viridis,
  clustering_method = "complete",#默认
  drop_levels = F,
  main = paste("Pancreatic cancer", " dependency matrix", sep = "")
)
dev.off()

# data.table::cbindlist
# ewid.data_gene.Sanger_bcscore_1=ewid.data_gene.Sanger_bcscore [names(ewid.data_gene.Sanger_bcscore)%in%screen_information_pancreas$EWHA_ID]
# ewid.data_gene.DepMap_20Q2_bcscore_1=ewid.data_gene.DepMap_20Q2_bcscore [names(ewid.data_gene.DepMap_20Q2_bcscore)%in%screen_information_pancreas$EWHA_ID]
# ewid.data_gene.GeCKO_bcscore_1=ewid.data_gene.GeCKO_bcscore [names(ewid.data_gene.GeCKO_bcscore)%in%screen_information_pancreas$EWHA_ID]
# ewid.data_gene.PubMed_1=ewid.data_gene.PubMed [names(ewid.data_gene.PubMed)%in%screen_information_pancreas$EWHA_ID]
ewid.data_gene_1=list(ewid.data_gene.Sanger_bcscore_1,ewid.data_gene.DepMap_20Q2_bcscore_1,ewid.data_gene.GeCKO_bcscore_1,ewid.data_gene.PubMed_1)
ewid.data_gene_1=unlist(ewid.data_gene_1,recursive =F)
ewid.data_gene_1[1:2]
ewid.data_gene_1=ewid.data_gene_1[names(ewid.data_gene_1)%in%names(resultTlist)]
load(file="C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/probe2gene_mrna_10.Rdata")
load(file="C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/probe2gene_mrna_25.Rdata")
ewid.data_gene_2=lapply(ewid.data_gene_1,function(x) subset(x,rownames(x)%in%probe2gene_mrna_10$symbol))
ewid.data_gene_2=lapply(ewid.data_gene_2,function(x) subset(x,select=3))
ewid.data_gene_2=lapply(ewid.data_gene_2,function(x) as.data.frame(t(x)))
ewid.data_gene_2=lapply(ewid.data_gene_2,function(x) subset(x,select=probe2gene_mrna_10$symbol))
ewid.data_gene_3=data.table::rbindlist(ewid.data_gene_2)
ewid.data_gene_3=as.data.frame(ewid.data_gene_3)
rownames(ewid.data_gene_3)=names(ewid.data_gene_2)
ewid.data_gene_3=t(ewid.data_gene_3)
ewid.data_gene_3=as.data.frame(ewid.data_gene_3)
writexl::write_xlsx(ewid.data_gene_3,"ewid.data_gene_3.xlsx")
##### 手动检查
gnenNA=colnames(ewid.data_gene_3)[which(apply(ewid.data_gene_3,2,function(x) all(is.na(x))))]
ewid.data_gene_3<-ewid.data_gene_3[,-which(apply(ewid.data_gene_3,2,function(x) all(is.na(x))))]
ewid.data_gene_3$median=apply(ewid.data_gene_3[,1:(nrow(ewid.data_gene_3))],1,function(x) median(x,na.rm=T))
ewid.data_gene_3=as.data.frame(ewid.data_gene_3)
ewid.data_gene_3=ewid.data_gene_3[order(ewid.data_gene_3$median,decreasing = F),]
ewid.data_gene_3$median=NULL
ewid.data_gene_4=as.data.frame(t(ewid.data_gene_3))
ewid.data_gene_4$EWHA_ID=rownames(ewid.data_gene_4)

screen_information_pancreas=as.data.frame(screen_information_pancreas)
library_information=as.data.frame(library_information)
cell_line_information_pancreas=as.data.frame(cell_line_information_pancreas)
ewid.data_gene_4=merge(ewid.data_gene_4,screen_information_pancreas,by="EWHA_ID")
ewid.data_gene_4=merge(ewid.data_gene_4,library_information,by="LIBRARY_ID")
ewid.data_gene_4=merge(ewid.data_gene_4,cell_line_information_pancreas,by="CELL_LINE")
nrow(ewid.data_gene_4)==ncol(ewid.data_gene_3)
ewid.data_gene_4=subset(ewid.data_gene_4,select = c(4:ncol(ewid.data_gene_4),1,2,3))
ewid.data_gene_4$source2=str_split(ewid.data_gene_4$SOURCE,"_",simplify = T)[,1]
ewid.data_gene_4=tidyr::unite(ewid.data_gene_4,"source_ccle_line",CELL_LINE,source2,sep="_",remove = F)
save(ewid.data_gene_4,ewid.data_gene_3,file = "ewid.data_gene_34.Rdata")
box=ewid.data_gene_4
library(tidyr)
box <- gather(box, key = "variable",value ="value", 1:(nrow(ewid.data_gene_3)) )
box=subset(box,!is.na(box$value))
library(ggpubr)

# google search : ggpubr boxplot add p-value
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
# ,outlier.shape = NA,order=

p <- ggboxplot(box, x = "variable", y = "value",  fill = "variable",outlier.shape = NA,#fill="variable",
               font.label = list(size = 0.4 )#,palette=c(1,2,3),order=c("No changes","Shallow Deletion","Deep Deletion")#"#00AFBB","#E7B800",
)+
  #add = "jitter"
  # scale_color_manual(values=c(1,2,3,4))  +#"#E69F00", "#56B4E9"
  coord_cartesian(ylim = c(-7.5, 2.5))+
  # theme(legend.position="left")+
  labs(title="",x="Batch corr", y = "Expression")+#rotate()
  rotate_x_text(90)
# ggpar(p,legend = "none")
#+#tag="ENSG00000000003"+
# ??geom_dotplot

p+geom_point(aes( color=source2),shape=1,binaxis='y', stackdir='center', stackratio=0.1, dotsize=0.5,show.legend = T)+  guides(color = guide_legend(title = 'Source',nrow = 1, byrow = TRUE,legend.position = "top", direction = "vertical"),fill="none")+
  theme(legend.position = c(0.8,0.1), legend.background = element_blank())

ggsave('gene_up_all_GO_dotplot.pdf',width = 9.6,height=6,dpi = 600)

# legend.position="top",
#  Add p-value
p + stat_compare_means(,method = "kruskal.test", label.y = 20 )#aes(group=cnv3)
# stat_compare_means(aes(group=cnv3),method = "kruskal.test", label.y = as.numeric(seq(14.7,200,by=0.7)))#+ #备注3

# ,comparisons =combn(1:5, 2, FUN = list)
#   stat_compare_means(comparisons = list(c("No changes","Shallow Deletion"),c("No changes","Deep Deletion"),c("Shallow Deletion","Deep Deletion")),method = "wilcox.test",label.y = as.numeric(seq(100,850,by=20)))+ #备注3
#   stat_summary(aes(x = cnv3, y = number1),fun= "mean", geom = "point",position = position_dodge(0.75), shape = 23, size = 1, fill = "pink")
# 
}

######

facp=fac
rownames(annotation_row) = paste("Gene", 1:20, sep = "")


annotation_row = data.frame(
  GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
)
rownames(annotation_row) = paste("Gene", 1:20, sep = "")


annotation_col = data.frame(
  "cell"=fac)
rownames(annotation_col) = colnames(dat)
ann_colors = list(
  "Studies" = c("Our data" = "#FF0000", "GSE97003" = "#0000FF", "CCLE" = "#006400")
)

#mhighall_over

library(cowplot)
library(pheatmap)
library(survival)
library(survminer)
library(data.table)
library(pheatmap)
library(gplots)
library(ggplot2)
library(pheatmap)
library(viridis)
library(edgeR)
exp_dat=ui[rownames(ui)%in%probe2gene_mrna_10$ENSEMBL,]

tmp=t(exp_dat)
tmp[tmp > 2] = 2
tmp[tmp < -2] = -2
tmp=t(tmp)
p =pheatmap(
  mat               = tmp,#t(cell_tumor_matrix)
  # color             = inferno(50),
  # border_color      = NA,
  show_colnames     = T,
  show_rownames     = T,
  fontsize          = 8,
  # annotation_col = annotation_col,
  annotation_legend = T,
  # annotation_colors = ann_colors,
  clustering_method = "complete",#默认
  drop_levels = F,
  main = paste("Pancreatic cancer", " dependency matrix", sep = "")
)
dev.off()
facp=as.data.frame(fac)
library(dplyr)
library(grid)
facp$colors = dplyr::case_when(fac=="Our data" ~ "#FF0000", fac =="GSE97003" ~ "#0000FF",fac =="CCLE" ~ "#006400")
cols = facp[order(match(colnames(dat), p$gtable$grobs[[5]]$label)),]$colors  #Assuming row labels are in grob 5
p$gtable$grobs[[5]]$gp = gpar(col = cols, fontface = "bold",fontsize= 8)
pdf(file='Irrr.pdf',height = 6,width = 12,onefile = T)
p
dev.off()
######################
new_dat=v0d
exp_dat=new_dat[names(sort(fp )),df$SYMBOL]#c('KIF15','FEN1','ZFP69B','SP6','SPARC'
tmp=t(exp_dat)
tmp[tmp > 2] = 2
tmp[tmp < -2] = -2
plot.h=pheatmap(tmp,show_colnames = F,cluster_cols = T,onefile = T)#
plot.h=pheatmap(main = " ",tmp,show_colnames = F,cluster_cols = F,cluster_rows = F,col= mycolors)#,col= mycolors
plot_grid(plot.point, plot.sur, plot.h$gtable,
          labels = c("", "",""),
          axis = "tb",
          align = 'v',ncol = 1)
dev.off()
pdf("Construction_risk_score.pdf",width = 8,height=8)
plot_grid(plot.point, plot.sur, plot.h$gtable,
          labels = c("","",""),
          axis = "tb",
          align = 'v',ncol = 1)
dev.off()
