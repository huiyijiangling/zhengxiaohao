#zxh
#20201018
project <- 'TCGA-PAAD'
clinicaldir <- paste(project, 'Clinical', sep='/')#这里该有的全有了
# gdcClinicalDownload(project.id     = 'TCGA-PAAD', 
#                     write.manifest = FALSE,
#                     method         = 'gdc-client',
#                     directory      = clinicaldir)
#By specifying key.info=TRUE, only common clinical information will be organized and reported. Otherwise, all the clinical information from the XML files will be extracted.
clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = F)
colnames(clinicalDa)<-colnames(clinicalDa)%>% stringr::str_replace_all("\\s+|\\.","-")#\\s 这就是单一个空格,把原来的点和空格变成-
clinicalDa[1:6,5:10]
clinicalDa_PAAD=as.data.frame(t(clinicalDa))
write.csv(clinicalDa,file="PAAD_clinical_gdc.csv",quote = T)
table(clinicalDa_PAAD$histological_type)
clinicalDa_PAAD=filter(clinicalDa_PAAD,histological_type=='Pancreas-Adenocarcinoma Ductal Type')
# clinicalDa_PAAD <- clinicalDa_PAAD[clinicalDa_PAAD$histological_type %in% "Pancreas-Adenocarcinoma Ductal Type",]
write.csv(clinicalDa_PAAD,file="PAAD_ductal_clinical_gdc.csv",quote = T)
save(clinicalDa_PAAD,file = "PAAD_ductal_clinical_gdc.Rdata")

# project <- 'TCGA-PAAD'
clinicaldir2 <- paste(project, 'bcr_biotab', sep='/')# 不能这么合并,剩下这七个需要分开合并
clinicalDa2 <- gdcClinicalMerge(path = clinicaldir2, key.info = F)
clinicalDa2[1:6,5:10]
write.csv(clinicalDa,file="PAAD_clinical7.csv",quote = T)

gsub(colnames(clinicalDa),)