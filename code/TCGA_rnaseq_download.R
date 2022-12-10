# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
gdcdata=function(i){
  library(TCGAbiolinks)
  projects <- getGDCprojects()
  library(dplyr)
  projects <- projects %>% 
    as.data.frame() %>% 
    select(project_id,tumor) %>% 
    filter(grepl(pattern="TCGA",project_id))
  ## 0.运行信息
  print(paste0("Di=1i=ownloading number ",i,",project name: ",projects$project_id[i]))
  ## 1.查询信息
  query.exp = GDCquery(project = projects$project_id[i], 
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")
  ## 2.正式下载
  GDCdownload(query.exp,
              method = "api", 
              files.per.chunk = 10)
  ## 3.多个数据合并
  pre.exp = GDCprepare(query = query.exp)
  ## 4.提取表达量数据
  library(SummarizedExperiment)
  countsdata = SummarizedExperiment::assay(pre.exp,1)
  tpmdata=SummarizedExperiment::assay(pre.exp,4)
  fpkmdata=SummarizedExperiment::assay(pre.exp,5)
  fpkmuqdata=SummarizedExperiment::assay(pre.exp,6)
  
  gene_id=data.frame(id=rowData(pre.exp)@listData[["gene_id"]], gene_name= rowData(pre.exp)@listData[["gene_name"]],gene_type=rowData(pre.exp)@listData[["gene_type"]])
  counts=cbind(gene_id,countsdata)
  tpm=cbind(gene_id,tpmdata)
  fpkm=cbind(gene_id,fpkmdata)
  fpkmuq=cbind(gene_id,fpkmuqdata)
  #临床信息
  clinical <- GDCquery_clinic(project = projects$project_id[i], type = "clinical")
  ## 5.保存数据
  filename1 = paste0("result/",projects$project_id[i],"-counts.Rdata")
  filename2 = paste0("result/",projects$project_id[i],"-tpm.Rdata")
  filename3 = paste0("result/",projects$project_id[i],"-fpkm.Rdata")
  filename4 = paste0("result/",projects$project_id[i],"-fpkmuq.Rdata")
  filename5 = paste0("result/",projects$project_id[i],"-clinical.Rdata")
  save(counts,file=filename1) 
  save(tpm,file=filename2)
  save(fpkm,file=filename3) 
  save(fpkmuq,file=filename4) 
  save(clinical,file=filename5) 
}
dir.create("result")
for (i in 1:33) {
  gdcdata(i)
}
##########################################
#cptac download
gdcdata=function(i){
  library(TCGAbiolinks)
  projects <- getGDCprojects()
  library(dplyr)
  projects <- projects %>% 
    as.data.frame() %>% 
    select(project_id,tumor) %>% 
    filter(grepl(pattern="CPTAC",project_id))
  ## 0.运行信息
  print(paste0("Di=1i=ownloading number ",i,",project name: ",projects$project_id[i]))
  ## 1.查询信息
  query.exp = GDCquery(project = projects$project_id[i], 
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       data.format="TSV",
                       workflow.type = "STAR - Counts")
  
  ## 2.正式下载
  GDCdownload(query.exp,
              method = "api", 
              files.per.chunk = 10)
  ## 3.多个数据合并
  pre.exp = GDCprepare(query = query.exp)
  ## 4.提取表达量数据
  library(SummarizedExperiment)
  countsdata = SummarizedExperiment::assay(pre.exp,1)
  tpmdata=SummarizedExperiment::assay(pre.exp,4)
  fpkmdata=SummarizedExperiment::assay(pre.exp,5)
  fpkmuqdata=SummarizedExperiment::assay(pre.exp,6)
  
  gene_id=data.frame(id=rowData(pre.exp)@listData[["gene_id"]], gene_name= rowData(pre.exp)@listData[["gene_name"]],gene_type=rowData(pre.exp)@listData[["gene_type"]])
  counts=cbind(gene_id,countsdata)
  tpm=cbind(gene_id,tpmdata)
  fpkm=cbind(gene_id,fpkmdata)
  fpkmuq=cbind(gene_id,fpkmuqdata)
  #临床信息
  clinical <- GDCquery_clinic(project = projects$project_id[i], type = "clinical")
  ## 5.保存数据
  filename1 = paste0("result/",projects$project_id[i],"-counts.Rdata")
  filename2 = paste0("result/",projects$project_id[i],"-tpm.Rdata")
  filename3 = paste0("result/",projects$project_id[i],"-fpkm.Rdata")
  filename4 = paste0("result/",projects$project_id[i],"-fpkmuq.Rdata")
  filename5 = paste0("result/",projects$project_id[i],"-clinical.Rdata")
  save(counts,file=filename1) 
  save(tpm,file=filename2)
  save(fpkm,file=filename3) 
  save(fpkmuq,file=filename4) 
  save(clinical,file=filename5) 
}
# cptac3
gdcdata=function(i){
  library(TCGAbiolinks)
  projects <- getGDCprojects()
  library(dplyr)
  projects <- projects %>% 
    as.data.frame() %>% 
    select(project_id,tumor) %>% 
    filter(grepl(pattern="CPTAC",project_id))
  ## 0.运行信息
  print(paste0("Di=1i=ownloading number ",i,",project name: ",projects$project_id[i]))
  ## 1.查询信息
  query.exp = GDCquery(project = projects$project_id[i], 
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       data.format="TSV",
                       workflow.type = "STAR - Counts")

  clinical_cptac3supp=query.exp$results[[1]]
  save(fefefe,file = "cptacfefefe.Rdata")#cptac3 无语了
  table(duplicated(query.exp$results[[1]]$sample.submitter_id))
  table(duplicated(query.exp$results[[1]]$cases.submitter_id))
  table(duplicated(query.exp$results[[1]]$cases))
  ## 2.正式下载
  GDCdownload(query.exp,
              method = "api", 
              files.per.chunk = 10)
# if(F){ 
#   
#   load(r"(C:\Users\zxh\Desktop\杜永星老师 提交\fto\cptac_data3.Rdata)")
#   query.exp$results[[1]]=subset(query.exp$results[[1]],!grepl(";",query.exp$results[[1]]$cases.submitter_id))#cptac3 only
#   query.exp$results[[1]]=subset(query.exp$results[[1]],!grepl(";",query.exp$results[[1]]$sample.submitter_id))#cptac3 only
#   query.exp$results[[1]]$cases.submitter_id=query.exp$results[[1]]$cases
#   query.exp$results[[1]]$sample.submitter_id=query.exp$results[[1]]$cases
# }
  ## 3.多个数据合并
  pre.exp = GDCprepare(query = query.exp,summarizedExperiment=F)
  pre.exp123 = pre.exp[,1:3]
  pre.exp4 = pre.exp[,-c(1:3)]
  unstranded
  ## 4.提取表达量数据
  library(SummarizedExperiment)
  
  countsdata = pre.exp4[,1:1883]
  tpmdata=subset(pre.exp4,select=grepl("tpm_unstranded",colnames(pre.exp4)))
  fpkmdata=subset(pre.exp4,select=grepl("fpkm_unstranded",colnames(pre.exp4)))
  fpkmuqdata=subset(pre.exp4,select=grepl("fpkm_uq_unstranded",colnames(pre.exp4)))
  
  gene_id=data.frame(pre.exp123)
  counts=cbind(gene_id,countsdata)
  tpm=cbind(gene_id,tpmdata)
  fpkm=cbind(gene_id,fpkmdata)
  fpkmuq=cbind(gene_id,fpkmuqdata)
  #临床信息
  clinical <- GDCquery_clinic(project = projects$project_id[i], type = "clinical")
  ## 5.保存数据
  filename1 = paste0("result/",projects$project_id[i],"-counts.Rdata")
  filename2 = paste0("result/",projects$project_id[i],"-tpm.Rdata")
  filename3 = paste0("result/",projects$project_id[i],"-fpkm.Rdata")
  filename4 = paste0("result/",projects$project_id[i],"-fpkmuq.Rdata")
  filename5 = paste0("result/",projects$project_id[i],"-clinical.Rdata")
  save(counts,file=filename1) 
  save(tpm,file=filename2)
  save(fpkm,file=filename3) 
  save(fpkmuq,file=filename4)
  table()
  clinical_cptac3supp$sample_type=stringr::str_split(clinical_cptac3supp$sample_type,";",simplify = T)[,1]
  save(clinical,clinical_cptac3supp,file=filename5) 
}





for (i in 1:2) {
  gdcdata(i)
}

### 作者大大的解决方案

# query.exp <- GDCquery(
#   project = 'CPTAC-3',
#   legacy = F,
#   data.category = "Transcriptome Profiling",
#   data.type = 'Gene Expression Quantification',
#   workflow.type = "STAR - Counts"
# )
# # remove duplicated
# query.exp$results[[1]] <- query.exp$results[[1]][!duplicated(query.exp$results[[1]]$sample.submitter_id),]
# 
# GDCdownload(query.exp,files.per.chunk = 40)
# se <- GDCprepare(
#   query = query.exp,
#   save = F
# )




# jjj=list.files(r"(C:\Users\zxh\Desktop\R\meta collect\GDCdata\CPTAC-3\harmonized\Transcriptome_Profiling\Gene_Expression_Quantification)",pattern = "tsv",recursive = T,
#                full.names = F)
# 
# table(duplicated(stringr::str_split(jjj,"\\/",simplify = T)[,2]))