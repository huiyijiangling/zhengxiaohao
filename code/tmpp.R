aa=readr::read_tsv(r"(C:\Users\zxh\Downloads\table (6).tsv)")
bb=readr::read_tsv(r"(C:\Users\zxh\Downloads\table (7).tsv)")
aa=aa[aa$`q-Value`<0.05,]
bb=bb[bb$`q-Value`<0.05,]
cc=merge(aa,bb,by="Correlated Gene")
cc=subset(cc,sign(cc$`Spearman's Correlation.x`)==sign(cc$`Spearman's Correlation.y`))
write.csv(cc,file="common_fto_rna_protein.csv")
  


