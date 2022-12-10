#1 6个方便的
query <-GDCquery(project = "TCGA-STAD", 
                 data.category = "Clinical",
                 data.type = "Clinical Supplement", 
                 # workflow.type,
                 legacy = FALSE,
                 access = "open",
                 # platform,
                 # file.type,
                 # barcode, 
                 # experimental.strategy
                 data.format = "bcr biotab"
                 # sample.type = c("Primary solid Tumor"),
                 )

GDCdownload(query)
clinical.BCRtab.all <-GDCprepare(query)
names(clinical.BCRtab.all)
save(clinical.BCRtab.all,file="clinical.BCRtab.all.Rdata")
load(file="clinical.BCRtab.all.Rdata")

#2 看不懂共14
query <-GDCquery(project = "TCGA-STAD", 
                 data.category = "Clinical",
                 data.type = "Clinical Supplement", 
                 # workflow.type,
                 legacy = FALSE,
                 access = "open",
                 # platform,
                 # file.type,
                 # barcode, 
                 # experimental.strategy
                 data.format = "bcr omf xml"
                 # sample.type = c("Primary solid Tumor"),
)

GDCdownload(query,method = "api",files.per.chunk = 5)
clinical.info<-c("drug","follow_up","radiation","patient","stage_event","new_tumor_event","admin")
for(i in clinical.info){
  clinical <- GDCprepare_clinic(query, clinical.info = i)
  write.csv(clinical, file = paste0('TCGA-STAD_clinical_',i,'.csv'), row.names =F)
}

#3 完整443

query <-GDCquery(project = "TCGA-STAD", 
                 data.category = "Clinical",
                 data.type = "Clinical Supplement", 
                 # workflow.type,
                 legacy = FALSE,
                 access = "open",
                 # platform,
                 # file.type,
                 # barcode, 
                 # experimental.strategy
                 data.format = "bcr xml"
                 # sample.type = c("Primary solid Tumor"),
)

GDCdownload(query,method = "api",files.per.chunk = 5)
clinical.info<-c("drug","follow_up","radiation","patient","stage_event","new_tumor_event","admin")
for(i in clinical.info){
  clinical <- GDCprepare_clinic(query, clinical.info = i)
  write.csv(clinical, file = paste0('TCGA-STAD_clinical_',i,'.csv'), row.names =F)
}

###############################
bcr auxiliary xml
443

bcr ssf xml
443

bcr xml
443

bcr biotab
11
?GDCprepare_clinic
query <-GDCquery(project = "TCGA-STAD", 
                 data.category = "Clinical", 
                 file.type = "xml")
GDCdownload(query)
clinical <-GDCprepare_clinic(query, clinical.info = "patient")  #上表中黄色部分可设置
#循环输出所有的临床数据：
clinical.info<-c("drug","follow_up","radiation","patient","stage_event","new_tumor_event","admin")
for(i in clinical.info){
  clinical <- GDCprepare_clinic(query, clinical.info = i)
  write.csv(clinical, file = paste0('TCGA-STAD_clinical_',i,'.csv'), row.names =F)
}