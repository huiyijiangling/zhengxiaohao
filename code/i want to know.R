load(r"(C:\Users\zxh\Desktop\R\meta collect\ensg_symbol_relation.Rdata)")
load(r"(C:\Users\zxh\Desktop\R\meta collect\tcga_pcawg_rnaseq.Rdata)")
load(r"(C:\Users\zxh\Desktop\R\meta collect\tcga_pcawg_rnaseq_counts.Rdata)")
pcawg_clinical_CDR=pcawg_clinical_CDR[pcawg_clinical_CDR$sample%in%phe_cachexia$sample,]

delist=dePC
eee=kmTestFun_single(rownames(delist),TCGA_tpm_ensg[["PAAD"]],phe_cachexia,sep='best')
eee0=eee

eee=coxphTestFun_single(rownames(delist),TCGA_tpm_ensg[["PAAD"]],phe_cachexia)
eee1=eee

# eee=coxphTestFun_single(rownames(delist),rnaCountslog21_Pancreas,phe_cachexia)
# eee2=eee
# eee2=cbind(eee0,eee1,eee2)