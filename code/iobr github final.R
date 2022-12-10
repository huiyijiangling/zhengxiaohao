rm(list=ls())
options(stringsAsFactors = F)
gc()

if(F){
  gmt=list.files(r"(C:\Users\zxh\Desktop\R\MSigDB\split gsea 202203\symbols)",pattern = ".gmt",full.names = T)
  gmt_name=list.files(r"(C:\Users\zxh\Desktop\R\MSigDB\split gsea 202203\symbols)",pattern = ".gmt",full.names = F)
  gmt_name=stringr::str_split(gmt_name,".v7.5.1.symbols.gmt",simplify = T)[,1]
  format_msigdb=lapply(gmt,function(x) format_msigdb(gmt=x, ont = "term", gene = "gene"))
  names(format_msigdb)=gmt_name
  #刘明洋导出
  if(F){
    iobr_sig_collect_nogene33=unlist(iobr_sig_collect_nogene)
    iobr_sig_collect_nogene33=unlist(iobr_sig_collect_nogene33)
    View(iobr_sig_collect_nogene33)
    openxlsx::write.xlsx(iobr_sig_collect_nogene,file = "iobr_sig_collect_nogene.xlsx")
    
    
    format_msigdbee=unlist(format_msigdb,recursive = F)
    table(grepl("cachexia",names(format_msigdbee),ignore.case = T))
    
    cachexia=format_msigdbee[grepl("cachexia",names(format_msigdbee),ignore.case = T)]
    openxlsx::write.xlsx(cachexia,file = "cachexia.xlsx")
    muscle=format_msigdbee[grepl("muscle",names(format_msigdbee),ignore.case = T)]
    muscle=as.matrix(muscle)
    muscle=as.data.frame(muscle)
    muscle$V2=rownames(muscle)
    openxlsx::write.xlsx(muscle,file = "muscle.xlsx")
    
  }
  save(format_msigdb,file="format_msigdb.Rdata")
  # format_msigdb=unlist(format_msigdb,recursive = F)
  # format_msigdb=data.table::rbindlist(format_msigdb)
  # sig_group属于内置已经算好了
  # 实际上需要我们从文章中自己去算的
data(package="IOBR")
#  "cellmarkers","BRef", "TRef",  "ips_gene_set "
iobr_sig_collect_name=c( "go_bp", "go_cc", "go_mf", "hallmark", "kegg", "reactome", "signature_collection", "signature_metabolism", "signature_tme", "signature_tumor")
data_frame_list=lapply(iobr_sig_collect_name,function(x) eval(parse(text=x)))
lapply(data_frame_list,function(x) table(duplicated(names(x))))#最里层没有重复
names(data_frame_list)=iobr_sig_collect_name
data_frame_list_all=unlist(data_frame_list,recursive = F)
# table(duplicated(data_frame_list_all))
# duplicated(list(A=c(1,2,3)), list(A=c(3,2,1)))#F
# table(duplicated(names(data_frame_list_all)))#F
data_frame_list_all_distinct=lapply(data_frame_list_all, list)

data_frame_list_all_distinct=tibble(names=names(data_frame_list_all_distinct),genelist=data_frame_list_all_distinct)
data_frame_list_all_distinct$name1=stringr::str_split(data_frame_list_all_distinct$names,"\\.",simplify = T,n = 2)[,1]
data_frame_list_all_distinct$name2=stringr::str_split(data_frame_list_all_distinct$names,"\\.",simplify = T,n = 2)[,2]
data_frame_list_all_distinct_nodupname=data_frame_list_all_distinct[data_frame_list_all_distinct$name2%in%data_frame_list_all_distinct[duplicated(data_frame_list_all_distinct$name2),]$name2,]
table(data_frame_list_all_distinct_nodupname$name1)
data_frame_list_all_distinct_nodupname_col=subset(data_frame_list_all_distinct_nodupname,name1%in%"signature_collection")

data_frame_list_all=data_frame_list_all[!names(data_frame_list_all)%in%data_frame_list_all_distinct_nodupname_col$names]

# tibble_fff=tibble(A=c("A","A"),B=c(list(A=list(1,2,3)), list(A=list(3,2,1))))
# setequal(tibble_fff$B)
library(AnnoProbe)
source(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\updateName.R)")
library(doParallel)
probe2gene=unique(unlist(data_frame_list_all))
# data_frame_list_all_ensg=foreach(x=data_frame_list_all,.combine = "c")  %dopar% {
forthetree = alias2SymbolUsingNCBImodify(alias = probe2gene)
probe2gene = cbind(probe2gene, forthetree)
probe2gene = probe2gene[!is.na(probe2gene$ENSEMBL), ]#我们不研究没有在ENSEMBL里面没有
probe2gene = unique(probe2gene[, c("probe2gene", "ENSEMBL")])#bioc
colnames(probe2gene) = c("probe_id", "ENSEMBL")
cl <- makeCluster(30)#makeCluster(detectCores())
registerDoParallel(cl)
data_frame_list_all_ensg=foreach(x=data_frame_list_all)  %dopar%{unique(probe2gene[probe2gene$probe_id %in%x,"ENSEMBL"])}

stopCluster(cl)
#没有损失
names(data_frame_list_all_ensg)=names(data_frame_list_all)
data_frame_list_all_ensg=data_frame_list_all_ensg[!is.na(data_frame_list_all_ensg)]

save(data_frame_list_all_ensg,data_frame_list_all,file="data_frame_list_all_ensg.Rdata")
}
if(F){
# names(signature_metabolism)[!names(signature_metabolism)%in%names(signature_collection)]
load(file="data_frame_list_all_ensg.Rdata")
data_frame_list_all=data_frame_list_all[names(data_frame_list_all_ensg)]
data_frame_list_all_ensg_left=lapply(data_frame_list_all_ensg,function(x) Reduce(intersect,list(rownames(ui_clust),x)))
table(lengths(data_frame_list_all_ensg_left)/lengths(data_frame_list_all)>=0.8)

iobr_sig_collect=data_frame_list_all_ensg_left[lengths(data_frame_list_all_ensg_left)/lengths(data_frame_list_all)>=0.8]

sig_d_ssgsea =calculate_sig_score(
  pdata = NULL,
  eset = ui_clust,
  signature = iobr_sig_collect,
  method = "ssgsea",#pca,
  mini_gene_count = 3,
  column_of_sample = "ID",
  print_gene_propotion = T,
  adjust_eset = FALSE,
  parallel.size = 1L,
  print_filtered_signatures = T)

sig_d_pca =calculate_sig_score(
  pdata = NULL,
  eset = ui_clust,
  signature = iobr_sig_collect,
  method = "pca",#pca,
  mini_gene_count = 3,
  column_of_sample = "ID",
  print_gene_propotion = T,
  adjust_eset = FALSE,
  parallel.size = 1L,
  print_filtered_signatures = T)
iobr_sig_collect_nogene=unique(colnames(sig_d_pca)[-c(1:2)],colnames(sig_d_pca)[-c(1:2)])
iobr_sig_collect_nogene=data.frame(dataset=stringr::str_split(iobr_sig_collect_nogene,"\\.",simplify = T)[,1],iobr_sig_collect_nogene=iobr_sig_collect_nogene)
iobr_sig_collect_nogene=split(iobr_sig_collect_nogene$iobr_sig_collect_nogene,iobr_sig_collect_nogene$dataset)
save(sig_d_pca,sig_d_ssgsea,iobr_sig_collect_nogene,iobr_sig_collect,file = "sig_ui_clust.Rdata")
}
if(F){
  grepl("ENSGxxxxx",iobr_sig_collect,ignore.case = T)

  names(iobr_sig_collect)[grepl("mhc|hla",names(iobr_sig_collect),ignore.case = T)]
  names(iobr_sig_collect)[grepl("t_cell",names(iobr_sig_collect),ignore.case = T)]
  names(iobr_sig_collect)[grepl("caf|fibro",names(iobr_sig_collect),ignore.case = T)]
  names(iobr_sig_collect)[grepl("treg",names(iobr_sig_collect),ignore.case = T)]
  names(iobr_sig_collect)[grepl("macro",names(iobr_sig_collect),ignore.case = T)]
  names(iobr_sig_collect)[grepl("regulat",names(iobr_sig_collect),ignore.case = T)&grepl("t_cell",names(iobr_sig_collect),ignore.case = T)]
}
load(file = "sig_ui_clust.Rdata")
#只处理免疫相关的开关
if(F){
  library(IOBR)
  
  immportfilenames <- list.files(path ="C:/Users/zxh/Desktop/R/immport/split/",pattern = ".txt",full.name=T)
  immportfilenames_short <- list.files(path ="C:/Users/zxh/Desktop/R/immport/split/",pattern = ".txt",full.name=F)
  immportfilenames_short=stringr::str_split(immportfilenames_short,"\\.",simplify = T)[,1]
  immportfiles=list()
  for (i in 1:length(immportfilenames)){immportfiles[[i]]=read.table(immportfilenames[[i]],sep="\t",stringsAsFactors = F,fill = TRUE,encoding = "UTF-8",header = T)
  immportfiles[[i]]$Immport_class=immportfilenames_short[[i]]
  }
  ImmPort_files=data.table::rbindlist(immportfiles)
  ImmPort_files_only_symbol=ImmPort_files[,-5]
  ImmPort_files_only_symbol=dplyr::distinct(ImmPort_files_only_symbol)
  # 总的少一些不知道为啥
  ImmPort_files_summary=read.table(r"(C:\Users\zxh\Desktop\R\immport\summary\GeneList.txt)",sep="\t",stringsAsFactors = F,fill = TRUE,encoding = "UTF-8",header = T)
  ImmPort_files$Immport_class[!ImmPort_files$Symbol%in%ImmPort_files_summary$Symbol]
  table(ImmPort_files_summary$Category)
  length(unique(ImmPort_files_summary$Symbol))
  # table(unlist(iobr_sig_collect_nogene) %in%names(iobr_sig_collect))
  library(AnnoProbe)
  source(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\updateName.R)")
  probe2gene=unique(ImmPort_files$Symbol)
  forthetree = alias2SymbolUsingNCBImodify(alias = probe2gene)
  probe2gene = cbind(probe2gene, forthetree)
  probe2gene = probe2gene[!is.na(probe2gene$ENSEMBL), ]#我们不研究没有在ENSEMBL里面没有
  ImmPort_files_ensg=merge(ImmPort_files,probe2gene,by.x="Symbol",by.y="probe2gene")
iobr_sig_collect=iobr_sig_collect[unlist(lapply(iobr_sig_collect, function(x) any(x %in%ImmPort_files_ensg$ENSEMBL)))]
for (i in 1:length(iobr_sig_collect_nogene)) {
iobr_sig_collect_nogene[[i]]=iobr_sig_collect_nogene[[i]][iobr_sig_collect_nogene[[i]] %in%names(iobr_sig_collect)]
}
lapply(iobr_sig_collect_nogene,function(x) x %in%names(iobr_sig_collect))
}

iobr_sig_collect_nogene=lapply(iobr_sig_collect_nogene, function(x) stringr::str_split(x,"\\.",n = 2,simplify = T)[,2])
if(T){
iobr_sig_collect_nogene_single=iobr_sig_collect_nogene  
iobr_sig_collect_nogene_single[["kegg1"]]=iobr_sig_collect_nogene_single[["kegg"]][1:60]
iobr_sig_collect_nogene_single[["kegg2"]]=iobr_sig_collect_nogene_single[["kegg"]][61:119]
iobr_sig_collect_nogene_single[c(1,2,3,5,6)]=NULL

# names(iobr_sig_collect_nogene_single)=NULL
# iobr_sig_collect_nogene_single=unlist(iobr_sig_collect_nogene_single)
# 
# iobr_sig_collect_nogene_single=data.frame(V1=iobr_sig_collect_nogene_single,V2=iobr_sig_collect_nogene_single)
# iobr_sig_collect_nogene_single=split(iobr_sig_collect_nogene_single,iobr_sig_collect_nogene_single$V2)
# iobr_sig_collect_nogene_single=lapply(iobr_sig_collect_nogene_single, function(x) x[[1]])

}

c1=stringr::str_split(colnames(sig_d_pca),"\\.",n = 2,simplify = T)[,2]
c1[1]="Index"
c1[2]="ID"
colnames(sig_d_pca)=c1
c2=stringr::str_split(colnames(sig_d_ssgsea),"\\.",n = 2,simplify = T)[,2]
c2[1]="ID"
c2[2]="Index"
colnames(sig_d_ssgsea)=c2
table(duplicated(lapply(names(iobr_sig_collect), function(x) stringr::str_split(x,"\\.",n = 2,simplify = T)[,2])))
names(iobr_sig_collect)=lapply(names(iobr_sig_collect), function(x) stringr::str_split(x,"\\.",n = 2,simplify = T)[,2])
# https://iobr.github.io/IOBR/IOBR-VIGNETTE.html
pdata_group<-imvigor210_pdata[!imvigor210_pdata$BOR_binary=="NA",c("ID","BOR","BOR_binary")]

# data("imvigor210_pdata")
# data("imvigor210_eset")
# data("imvigor210_sig")
# 
# names(sig_group) %in% names(signature_collection)#F
# names(sig_group) 


#1
names(signature_collection)[!names(signature_collection)%in%unique(unlist(sig_group))]
#2
unique(unlist(sig_group))[!unique(unlist(sig_group))%in%names(signature_collection)]#少了
#3
data(sig_group)

sig_group$tumor_signature[1]
aaaaaa=unique(unlist(sig_group))
bbbbbb=unique(unlist(signature_collection))
names(signature_tme) %in% names(signature_collection)#T

#综上所述 signature_collection 这可以使用，但是离sig_group还少了点，补计算
# [1] "TMEscore_CIR"                                      
# [2] "IPS_IPS"          
# 没找到

unlist(iobr_sig_collect_nogene)%in%colnames(sig_d_ssgsea)

# pan_cluster3=read.csv(r"(C:\Users\zxh\Desktop\R\paad-tcga-gtex\crispr25 for final\crispr10.k=4.consensusClass.csv)",header =  F)
load(r"(C:\Users\zxh\Desktop\R\paad-tcga-gtex\dataset5 for clust.Rdata)")
ui_clust_da=as.data.frame(ui_clust)

if(F){
ui_clust_da=ui_clust_da[c("ENSG00000181847","ENSG00000135077","ENSG00000120217","ENSG00000197646","ENSG00000188389","ENSG00000089692","ENSG00000104760","ENSG00000107738","ENSG00000163599"),]
ui_clust_da=ui_clust_da[c("ENSG00000164530","ENSG00000143196","ENSG00000172061","ENSG00000143632","ENSG00000172061","ENSG00000197635"),]
}
# ,"ENSG00000104760","ENSG00000107738","ENSG00000163599"
# PI16
# Ensembl ID: ENSG00000164530
# DPT
# Ensembl ID: ENSG00000143196.4
# LRRC15
# Ensembl ID: ENSG00000172061.8
# ACTA1
# Ensembl ID: ENSG00000143632.14
# LRRC15
# Ensembl ID: ENSG00000172061.8
# DPP4
# Ensembl ID: ENSG00000197635.9

pan_cluster3=read.csv(r"(C:\Users\zxh\Desktop\R\paad-tcga-gtex\group_into_3.csv)",header =T)
pan_cluster3=tibble(pan_cluster3)
pan_cluster3$x=paste0("Cluster",pan_cluster3$x)
#try 这个函数真牛逼
#t
try({res1 <- iobr_cor_plot(
  pdata_group = pan_cluster3,
  id1                   = "X",
  feature_data          = sig_d_ssgsea,
  id2                   = "ID",
  target                = NULL,
  group                 = "x", #分类变量，如果是连续变量，group2或group3
  is_target_continuous  = FALSE,
  padj_cutoff = 0.05,
  index = 1,
  category = "signature",#'gene',属于'signature_collection' or 'signature_tme';;;#'signature',属于'sig_group'
  signature_group = iobr_sig_collect_nogene_single,
  #iobr_sig_collect_nogene[c(4)],#	'sig_group' , 'signature_collection' or 'signature_tme'
  ProjectID = "PAN",#这个我不确定
  palette_box = "set2",#"jco","paired1",#"nrc",
  cols_box = NULL,
  palette_corplot = "pheatmap",
  palette_heatmap = 2,#234都有
  feature_limit = 100000,#这个要比单个characeter多
  character_limit = 3000,#这个要大
  show_heatmap_col_name = FALSE,
  show_col = FALSE,
  show_plot = F,
  path = "d2",
  discrete_x = 20,
  discrete_width = 100,
  show_palettes = FALSE,
  fig.type = "pdf")},silent = T)

# pan_cluster37777777=pan_cluster3
# pan_cluster37777777[pan_cluster37777777$x=="Cluster2","x"]="Cluster3"
source(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\updateName.R)")
probe2gene=data.frame(probe_id=rownames(ui_clust_da),symbol=rownames(ui_clust_da))
forthetree=ensg2SymbolUsingNCBImodify(alias=rownames(ui_clust_da),gene.info.file = "Homo_sapiens.gene_info.gz")
probe2gene=cbind(probe2gene,forthetree)
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
probe2gene=unique(probe2gene[,c("probe_id","Symbol")])
# probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
colnames(probe2gene)=c("probe_id","symbol")
ui_clust_da <- filterEM(ui_clust_da,probe2gene)


library(doParallel)
probe2gene=unique(unlist(iobr_sig_collect))
names_iobr_sig_collect=names(iobr_sig_collect)
# data_frame_list_all_ensg=foreach(x=data_frame_list_all,.combine = "c")  %dopar% {
forthetree = ensg2SymbolUsingNCBImodify(alias = probe2gene)
probe2gene = cbind(probe2gene, forthetree)
probe2gene = probe2gene[!is.na(probe2gene$ENSEMBL), ]#我们不研究没有在ENSEMBL里面没有
probe2gene = unique(probe2gene[, c("probe2gene","Symbol")])#bioc
colnames(probe2gene) = c("probe_id", "symbol")
cl <- makeCluster(30)#makeCluster(detectCores())
registerDoParallel(cl)
iobr_sig_collect_symbol=foreach(x=iobr_sig_collect)  %dopar%{unique(probe2gene[probe2gene$probe_id %in%x,"symbol"])}
names(iobr_sig_collect_symbol)=names_iobr_sig_collect
iobr_sig_collect_symbol=iobr_sig_collect_symbol[!is.na(iobr_sig_collect_symbol)]
stopCluster(cl)

wwww=list(rownames(ui_clust_da))
names(wwww)="hotgene"


names(iobr_sig_collect)[grepl("t_cell",names(iobr_sig_collect),ignore.case = T)]

names(iobr_sig_collect)[grepl("t_cell",names(iobr_sig_collect),ignore.case = T)&(!grepl("nk",names(iobr_sig_collect),ignore.case = T))]
names(iobr_sig_collect)[grepl("caf|fibro",names(iobr_sig_collect),ignore.case = T)]

names(iobr_sig_collect)[grepl("macro",names(iobr_sig_collect),ignore.case = T)]

names(iobr_sig_collect)[grepl("treg",names(iobr_sig_collect),ignore.case = T)|(grepl("Regulatory",names(iobr_sig_collect),ignore.case = T)&grepl("t_cell",names(iobr_sig_collect),ignore.case = T))]

names(iobr_sig_collect)[grepl("b_cell",names(iobr_sig_collect),ignore.case = T)]




res2 <-iobr_cor_plot(
  pdata_group = pan_cluster3,
  id1 = "X",
  feature_data = ui_clust_da,#imvigor210_sig,
  # id2 = "ID",
  target = NULL,#连续变量
  group = "x", #分类变量，如果是连续变量，group2或group3
  is_target_continuous = F,#看情况
  padj_cutoff = 0.05,
  index = 1,
  category = "gene",#'gene',属于'signature_collection' or 'signature_tme';;;#'signature',属于'sig_group'
  signature_group = wwww,#iobr_sig_collect_symbol[grepl("FIBROBLAST|caf|fibro",names(iobr_sig_collect),ignore.case = T)],#wwww,#iobr_sig_collect_symbol[grepl("FIBROBLAST|caf|fibro",names(iobr_sig_collect),ignore.case = T)],
  
  #iobr_sig_collect[grepl("mhc",names(iobr_sig_collect),ignore.case = T)],#iobr_sig_collect,#	'sig_group' , 'signature_collection' or 'signature_tme'
  ProjectID = "PAN",#这个我不确定
  palette_box = "set2",#"jco","paired1",#"nrc",
  cols_box = NULL,
  palette_corplot = "pheatmap",
  palette_heatmap = 2,#234都有
  feature_limit = 260,#26
  character_limit = 600,#60
  show_heatmap_col_name = FALSE,
  show_col = FALSE,
  show_plot = T,
  path = "caf2",#"hotgene"
  discrete_x = 20,
  discrete_width = 20,
  show_palettes = FALSE,
  fig.type = "pdf")



# res<-
  iobr_cor_plot(pdata_group           = pan_cluster3,
                   id1                   = "V1",
                   feature_data          = sig_d_ssgsea,#ui_clust_da,
                   id2                   = "ID",
                   target                = NULL,
                   group                 = "V2",
                   is_target_continuous  = FALSE,
                   padj_cutoff           = 1,
                   index                 = 1,
                   category              = "signature",
                   signature_group       = iobr_sig_collect_nogene,#iobr_sig_collect,    
                   ProjectID             = "IMvigor210",
                   palette_box           = "set2",
                   palette_corplot       = "pheatmap",
                   palette_heatmap       = 2,
                   feature_limit         = 2600000,
                   character_limit       = 3000000,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = TRUE,
                   path                  = "2-BOR-relevant-genes")

  
  
  
  
  
  
  
  #####
  
  
  # commmon 通路 ras 上下
  load(file = "sig_ui_clust.Rdata")
  
  rasupdown=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\paad-tcga-gtex\1mapk2pi3k2kt3all.xlsx)",sheet = 4)
  rasupdown$GENEID=as.numeric(rasupdown$GENEID)
  library(AnnoProbe)
  library(IOBR) #加载IOBR
  source("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/updateName.R")
  
  dataset=rasupdown
  probe2gene=data.frame(probe_id=dataset,symbol=dataset)
  forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol,gene.info.file = "Homo_sapiens.gene_info.gz")
  probe2gene=cbind(probe2gene,forthetree)
  
  writexl::write_xlsx(EEEE,"clust gene list.xlsx")
  
  
  
  
  
  
  #########################################immport
  read