aa=openxlsx::read.xlsx(r"(C:\Users\zxh\Desktop\R\meta collect\xena\pcawg_specimen_histology_August2016_v9.xlsx)")
aa$level_of_cellularity=NULL
aa$percentage_cellularity=NULL
aa$tumour_grade=NULL
aa$tumour_stage=NULL
aa$tumour_histological_code=NULL
aa$icgc_specimen_id=NULL
aa=unique(aa)
apply(aa, 2, table)
openxlsx::write.xlsx(aa,file="historead.xlsx")
