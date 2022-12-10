#补齐一下
if(F){
  library(AnnoProbe)
  library(IOBR) #加载IOBR
  my_signature<-readxl::read_excel("C:/Users/zxh/Desktop/subtype/subtype clsut gene.xlsx")
  # my_signature[1:10,1:5]
  my_signature=my_signature[,c(1,3,5)]
  my_signature<-format_signatures(my_signature)
  EEEE <- list()
  source("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/updateName.R")
  for(i in 1:length(my_signature)){
  dataset=na.omit(my_signature[[i]])
  probe2gene=data.frame(probe_id=dataset,symbol=dataset)
  forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol,gene.info.file = "Homo_sapiens.gene_info.gz")
  probe2gene=cbind(probe2gene,forthetree)
  EEEE[[i]]=probe2gene}
  writexl::write_xlsx(EEEE,"clust gene list.xlsx")
}
if(T){
  # NA补齐的阵列
  # a b
  # hh hh
  # NA PDL1
  my_signature<-readxl::read_excel("C:/Users/zxh/Desktop/subtype/subtype clsut gene.xlsx") 
  my_signature[1:10,1:5]
  my_signature=my_signature[,c(2,5,9)]
  my_signature<-format_signatures(my_signature)
  my_signature[1:4]
  

}
load("C:/Users/zxh/Desktop/R/paad-tcga-gtex/dataset5 for clust.Rdata")

# rnatpm_Pancreas  
library(AnnoProbe)
dataset=ui_clust[,merrrr$sample]

probe2gene=data.frame(probe_id=rownames(dataset),symbol=rownames(dataset))
source("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/updateName.R")
forthetree=ensg2SymbolUsingNCBImodify(alias=probe2gene$symbol,gene.info.file = "Homo_sapiens.gene_info.gz")
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
if(F){
  forthetree=updateName(probe2gene$symbol)
  table(forthetree$gene_biotype)
  probe2gene=merge(forthetree,probe2gene,by.x="ALIAS",by.y="symbol",all.y=T)
}
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
probe2gene$GeneID=as.numeric(probe2gene$GeneID)
if(F){
  #ras
  probe2gene_ras3=probe2gene[probe2gene$GeneID %in%rasupdown$GENEID,]
  probe2gene_ras3=unique(probe2gene_ras3[,c("probe_id","Symbol")])
  # probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
  colnames(probe2gene_ras3)=c("probe_id","symbol")
  probes_expr=dataset
  genes_expr_ras3 <- filterEM(probes_expr,probe2gene_ras3)
}
probe2gene=unique(probe2gene[,c("probe_id","Symbol")])
# probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
colnames(probe2gene)=c("probe_id","symbol")
probes_expr=dataset
genes_expr <- filterEM(probes_expr,probe2gene)
# tide=genes_expr
# tide=t(scale(t(log2(tide+1)),center=T,scale = F))
# psych::describe(tide)
# write.table(tide,gzfile("newfile.txt.gz"),sep = "\t",quote = F)
# eset_stad 替换 genes_expr


library(IOBR) #加载IOBR
library(tidyr)
#TME associated signatures 
names(signature_tme)
View(signature_tme)
#Metabolism related signatures
names(signature_metabolism)

names(signature_tumor)


names(signature_collection)

signature_collection_citation
View(signature_collection_citation)
# data("eset_stad")
eset_stad[1:5, 1:5]

sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = genes_expr,
                             signature       = signature_collection,
                             method          = "pca",
                             mini_gene_count = 2)

tme_deconvolution_methods
#>         MCPcounter               EPIC              xCell          CIBERSORT 
#>       "mcpcounter"             "epic"            "xcell"        "cibersort" 
#> CIBERSORT Absolute                IPS           ESTIMATE                SVM 
#>    "cibersort_abs"              "ips"         "estimate"              "svm" 
#>               lsei              TIMER          quanTIseq 
#>             "lsei"            "timer"        "quantiseq"
# Return available parameter options of TME deconvolution.

help(deconvo_tme)
cibersort<-deconvo_tme(eset = genes_expr, method = "cibersort", arrays = FALSE, perm = 1000 )

#> 
#> >>> Running CIBERSORT
head(cibersort)
#> # A tibble: 6 x 26
#>   ID    B_cells_naive_C… B_cells_memory_… Plasma_cells_CI… T_cells_CD8_CIB…
#>   <chr>            <dbl>            <dbl>            <dbl>            <dbl>
#> 1 TCGA…           0.0323                0          0                 0.193 
#> 2 TCGA…           0.0866                0          0                 0.0872
#> 3 TCGA…           0.0474                0          0.00644           0.0286
#> 4 TCGA…           0.0125                0          0.00257           0.224 
#> 5 TCGA…           0.0544                0          0.00923           0.0936
#> 6 TCGA…           0.0246                0          0.00162           0.124 
#> # … with 21 more variables: T_cells_CD4_naive_CIBERSORT <dbl>,
#> #   T_cells_CD4_memory_resting_CIBERSORT <dbl>,
#> #   T_cells_CD4_memory_activated_CIBERSORT <dbl>,
#> #   T_cells_follicular_helper_CIBERSORT <dbl>,
#> #   `T_cells_regulatory_(Tregs)_CIBERSORT` <dbl>,
#> #   T_cells_gamma_delta_CIBERSORT <dbl>, NK_cells_resting_CIBERSORT <dbl>,
#> #   NK_cells_activated_CIBERSORT <dbl>, Monocytes_CIBERSORT <dbl>,
#> #   Macrophages_M0_CIBERSORT <dbl>, Macrophages_M1_CIBERSORT <dbl>,
#> #   Macrophages_M2_CIBERSORT <dbl>, Dendritic_cells_resting_CIBERSORT <dbl>,
#> #   Dendritic_cells_activated_CIBERSORT <dbl>,
#> #   Mast_cells_resting_CIBERSORT <dbl>, Mast_cells_activated_CIBERSORT <dbl>,
#> #   Eosinophils_CIBERSORT <dbl>, Neutrophils_CIBERSORT <dbl>,
#> #   `P-value_CIBERSORT` <dbl>, Correlation_CIBERSORT <dbl>,
#> #   RMSE_CIBERSORT <dbl>
pdf("/cell_bar_plot.pdf",width=8,height=8)
res<-cell_bar_plot(input = cibersort, title = "CIBERSORT Cell Fraction",coord_filp = F)
dev.off()
#> There are seven categories you can choose: box, continue2, continue, random, heatmap, heatmap3, tidyheatmap


epic<-deconvo_tme(eset = genes_expr, method = "epic", arrays = FALSE)#有报错
mcp<-deconvo_tme(eset = genes_expr, method = "mcpcounter")
xcell<-deconvo_tme(eset = genes_expr, method = "xcell",arrays = FALSE)
estimate<-deconvo_tme(eset = genes_expr, method = "estimate")
timer<-deconvo_tme(eset = genes_expr, method = "timer", group_list = rep("paad",dim(genes_expr)[2]))
quantiseq<-deconvo_tme(eset = genes_expr, tumor = TRUE, arrays = FALSE, scale_mrna = TRUE, method = "quantiseq")
ips<-deconvo_tme(eset = genes_expr, method = "ips", plot= FALSE)


######

library("IOBR")
library("tidyHeatmap")
data("imvigor210_sig")
imvigor210_sig[1:5, 1:10]
#> # A tibble: 5 x 10
#>   ID    B_cells_naive_C… B_cells_memory_… Plasma_cells_CI… T_cells_CD8_CIB…
#>   <chr>            <dbl>            <dbl>            <dbl>            <dbl>
#> 1 SAMf…           0.0290          0.0457           0                      0
#> 2 SAM6…           0.0812          0.00127          0                      0
#> 3 SAMc…           0.0124          0                0.00127                0
#> 4 SAM8…           0               0                0.00119                0
#> 5 SAMf…           0               0                0.00949                0
#> # … with 5 more variables: T_cells_CD4_naive_CIBERSORT <dbl>,
#> #   T_cells_CD4_memory_resting_CIBERSORT <dbl>,
#> #   T_cells_CD4_memory_activated_CIBERSORT <dbl>,
#> #   T_cells_follicular_helper_CIBERSORT <dbl>,
#> #   `T_cells_regulatory_(Tregs)_CIBERSORT` <dbl>

# Check therapeutic response and survival outcome in the phenotype data.
data("imvigor210_pdata")
imvigor210_pdata[1:5, 1:5]
#> # A tibble: 5 x 5
#>   ID              BOR   BOR_binary OS_days            OS_status
#>   <chr>           <chr> <chr>      <chr>              <chr>    
#> 1 SAM00b9e5c52da9 NA    NA         57.166324439999997 1        
#> 2 SAM0257bbbbd388 SD    NR         469.15811100000002 1        
#> 3 SAM025b45c27e05 PD    NR         263.16221766000001 1        
#> 4 SAM032c642382a7 PD    NR         74.907597539999998 1        
#> 5 SAM04c589eb3fb3 NA    NA         20.698151939999999 0

View(sig_group)
# The signature group list. 39个分组的名字
names(sig_group)
#>  [1] "tumor_signature"        "EMT"                    "io_biomarkers"         
#>  [4] "immu_microenvironment"  "immu_suppression"       "immu_exclusion"        
#>  [7] "immu_exhaustion"        "TCR_BCR"                "tme_signatures1"       
#> [10] "tme_signatures2"        "Bcells"                 "Tcells"                
#> [13] "DCs"                    "Macrophages"            "Neutrophils"           
#> [16] "Monocytes"              "CAFs"                   "NK"                    
#> [19] "tme_cell_types"         "CIBERSORT"              "MCPcounter"            
#> [22] "EPIC"                   "xCell"                  "quanTIseq"             
#> [25] "ESTIMATE"               "IPS"                    "TIMER"                 
#> [28] "fatty_acid_metabolism"  "hypoxia_signature"      "cholesterol_metabolism"
#> [31] "Metabolism"             "hallmark"               "hallmark1"             
#> [34] "hallmark2"              "hallmark3"              "Rooney_et_al"          
#> [37] "Bindea_et_al"           "Li_et_al"               "Peng_et_al"

# The tumor-relevant signatures in first group.
names(sig_group)[1]
#> [1] "tumor_signature"
sig_group[[1]]
#>  [1] "CellCycle_Reg"                            
#>  [2] "Cell_cycle"                               
#>  [3] "DDR"                                      
#>  [4] "Mismatch_Repair"                          
#>  [5] "Histones"                                 
#>  [6] "Homologous_recombination"                 
#>  [7] "Nature_metabolism_Hypoxia"                
#>  [8] "Molecular_Cancer_m6A"                     
#>  [9] "MT_exosome"                               
#> [10] "Positive_regulation_of_exosomal_secretion"
#> [11] "Ferroptosis"                              
#> [12] "EV_Cell_2020"

# The signatures associated with immunotherapy biomarkers.
names(sig_group)[3]
#> [1] "io_biomarkers"
sig_group[[3]]
#>  [1] "TMEscore_CIR"                    "TMEscoreA_CIR"                  
#>  [3] "TMEscoreB_CIR"                   "T_cell_inflamed_GEP_Ayers_et_al"
#>  [5] "CD_8_T_effector"                 "IPS_IPS"                        
#>  [7] "Immune_Checkpoint"               "Exhausted_CD8_Danaher_et_al"    
#>  [9] "Pan_F_TBRs"                      "Mismatch_Repair"                
#> [11] "APM"

# The signatures of immunosuppression.
names(sig_group)[5]
#> [1] "immu_suppression"
sig_group[[5]]
#> [1] "Pan_F_TBRs"                          
#> [2] "Fibroblasts_MCPcounter"              
#> [3] "Immune_Checkpoint"                   
#> [4] "Exhausted_CD8_Danaher_et_al"         
#> [5] "MDSC_Wang_et_al"                     
#> [6] "Macrophages_M2_cibersort"            
#> [7] "Tregs_quantiseq"                     
#> [8] "T_cells_regulatory_(Tregs)_CIBERSORT"
View(imvigor210_sig)
pdata_group<-imvigor210_pdata[!imvigor210_pdata$BOR_binary=="NA",c("ID","BOR","BOR_binary")]
res<-iobr_cor_plot(pdata_group           = pdata_group, #患者的分组数据
                   id1                   = "ID", #pdata的患者ID
                   feature_data          = imvigor210_sig, #患者解析出来的数据
                   id2                   = "ID", #imvigor210_sig的患者ID
                   target                = NULL, #目标还可以是连续性变量，比如你鉴定了一个新发现的signature，想要去看看它和已知的signatures的关系
                   group                 = "BOR_binary", #需要比较的分组，目前只支持二分类，以后会完善多个分组
                   is_target_continuous  = FALSE, #如果楼上的target是连续性变量==TRUE
                   padj_cutoff           = 1, #是否过滤变量
                   index                 = 1, #用来命名文件
                   category              = "signature", #这个功能还可以输入基因表达矩阵，来查看与免疫治疗反应相关的基因==“gene”
                   signature_group       = sig_group,#[c(1,3,5)], #这里选择了sig_group的1,3,5个分组进行可视化分析
                   ProjectID             = "IMvigor210", #用于文件的命名，避免结果文件的覆盖
                   palette_box           = "paired1", #选择boxplot的颜色，我们提供了好几个，准备添加自己命名颜色的功能；
                   palette_corplot       = "pheatmap",#选择corplot的颜色，我们提供了好几个，准备添加自己命名颜色的功能；
                   palette_heatmap       = 2, #选择heatmap的颜色，我们提供了好几个，准备添加自己命名颜色的功能；
                   feature_limit         = 26, #箱图展示变量的最大值/2
                   character_limit       = 30, #变量字符的最大值，大于30的字符会被去掉，防止影响图的美观
                   show_heatmap_col_name = FALSE,#是否展示热图的每个观测的名字
                   show_col              = FALSE, #是否显示画板
                   show_plot             = TRUE, #是否将图打印出来
                   path                  = "1-BOR-relevant-signatures")  #结果储存的路径



View(sig_group)
#> [1] ">>>  Processing signature: tumor_signature"
head(res)

library("IOBR")
# Load the test data: the gene expression matrix of IMvigor210 cohort has been normalized using method `voom`.
imvigor210_eset[1:5,1:5]
#>              SAMf2ce197162ce SAM698d8d76b934 SAMc1b27bc16435 SAM85e41e7f33f9
#> LOC100093631       3.1009398       2.8207636       3.7058993      2.81012848
#> LOC100126784      -2.6237747      -4.2560520      -5.4104447     -2.07600356
#> LOC100128108      -1.5017841       0.5200520      -0.7665885     -2.07600356
#> LOC100128288      -0.3361981      -1.2204281      -1.9510131     -1.25886761
#> LOC100128361       0.2545468       0.2923847      -0.2009913     -0.02537748
#>              SAMf275eb859a39
#> LOC100093631       4.0102463
#> LOC100126784      -3.6376118
#> LOC100128108      -1.9495558
#> LOC100128288       0.3320146
#> LOC100128361       0.7698920
View(imvigor210_eset)
# Extract the significant genes as a phenotype of the patients based on the exploration of expression matrix.
pdata_group<-rownames_to_column(as.data.frame(t(imvigor210_eset)),var = "ID")

pdata_group<-as.data.frame(pdata_group[,c("ID","HCP5","LINC00657")])
head(pdata_group)
#>                ID       HCP5 LINC00657
#> 1 SAMf2ce197162ce -1.8532565  2.746143
#> 2 SAM698d8d76b934 -1.7199991  1.636339
#> 3 SAMc1b27bc16435 -2.6030898  3.863351
#> 4 SAM85e41e7f33f9 -0.5375836  2.989060
#> 5 SAMf275eb859a39  0.4685876  2.599731
#> 6 SAM7f0d9cc7f001 -1.4439383  2.029233

imvigor210_sig[1:5, 1:5]

View(res)
View(sig_group)
res<-iobr_cor_plot(pdata_group           = pdata_group, #基因信息，连续性
                   id1                   = "ID", #pdata 的患者id
                   feature_data          = imvigor210_sig, #肿瘤微环境信息
                   id2                   = "ID", #feature_data的患者id
                   target                = "HCP5", #目标基因的列名
                   group                 = "group3", #这里详细解释一下这个分组的意思：见下文
                   is_target_continuous  = TRUE, #目标是否为连续性变量；
                   padj_cutoff           = 1, #这里不设阈值，否则会过滤掉很多变量
                   index                 = 3, #为了命名，可以不设置
                   category              = "signature", #因为这里的目标是通过基因来找相关的signature和细胞，所以设置成“signature”
                   signature_group       = sig_group[1:2], #signature的分组信息,是一个list
                   ProjectID             = "IMvigor210", #为了命名储存路径，可以不设置
                   palette_box           = "set2", #以下几个参数请看前期的推文
                   palette_corplot       = "pheatmap",
                   palette_heatmap       = 2,
                   feature_limit         = 26,
                   character_limit       = 30,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = TRUE,
                   path                  = "HCP5-relevant-signatures")



library("IOBR")
# help("make_mut_matrix")

maf_file<-"TCGA.STAD.mutect.c06465a3-50e7-46f7-b2dd-7bd654ca206b.DR-10.0.somatic.maf"
mut_list<-make_mut_matrix(maf      = maf_file,
                          isTCGA   = T, 
                          category = "multi")
#> -Reading
#> -Validating
#> --Removed 2 duplicated variants
#> -Silent variants: 70966 
#> -Summarizing
#> --Possible FLAGS among top ten genes:
#>   TTN
#>   MUC16
#>   SYNE1
#>   FLG
#> -Processing clinical data
#> --Missing clinical data
#> -Finished in 14.5s elapsed (25.5s cpu) 
#>        Frame_Shift_Del        Frame_Shift_Ins           In_Frame_Del 
#>                  18418                   4461                    692 
#>           In_Frame_Ins      Missense_Mutation      Nonsense_Mutation 
#>                    268                 109668                   6011 
#>       Nonstop_Mutation            Splice_Site Translation_Start_Site 
#>                    107                   2445                    106 
#>    DEL    INS    SNP 
#>  19387   4900 117889
#> Aggregation function missing: defaulting to length
#> Aggregation function missing: defaulting to length
#> Aggregation function missing: defaulting to length
#> Aggregation function missing: defaulting to length
# If "multi" is set in above "category" parameter, four data frames will be returned, which evaluate all the mutations of every gene or estimate only SNP, indel and frameshift as follow:

library("UCSCXenaTools")
var_stad<-XenaGenerate(subset = XenaCohorts =="GDC TCGA Stomach Cancer (STAD)") %>% 
  XenaFilter(filterDatasets    = "TCGA-STAD.mutect2_snv.tsv") %>% 
  XenaQuery() %>%
  XenaDownload() %>% 
  XenaPrepare()
#> This will check url status, please be patient.
#> All downloaded files will under directory /tmp/Rtmp0RFeyC.
#> The 'trans_slash' option is FALSE, keep same directory structure as Xena.
#> Creating directories for datasets...
#> Downloading TCGA-STAD.mutect2_snv.tsv.gz
head(var_stad)
#> # A tibble: 6 x 11
#>   Sample_ID gene  chrom  start    end ref   alt   Amino_Acid_Chan… effect filter
#>   <chr>     <chr> <chr>  <dbl>  <dbl> <chr> <chr> <chr>            <chr>  <chr> 
#> 1 TCGA-CD-… C1or… chr1  2.19e6 2.19e6 G     -     p.P72Rfs*87      frame… PASS  
#> 2 TCGA-CD-… ERRF… chr1  8.01e6 8.01e6 C     T     p.P327P          synon… panel…
#> 3 TCGA-CD-… CLCN6 chr1  1.18e7 1.18e7 G     A     p.S486N          misse… PASS  
#> 4 TCGA-CD-… PRAM… chr1  1.28e7 1.28e7 G     A     p.G341R          misse… panel…
#> 5 TCGA-CD-… PRAM… chr1  1.32e7 1.32e7 G     T     p.P148H          misse… PASS  
#> 6 TCGA-CD-… CELA… chr1  1.55e7 1.55e7 G     A     p.G59R           misse… PASS  
#> # … with 1 more variable: dna_vaf <dbl>
# Then function `make_mut_matrix` can be used to transform data frame into mutation matrix
mut_list2<-make_mut_matrix(mut_data               = var_stad,
                           category               = "multi",
                           Tumor_Sample_Barcode   = "Sample_ID",
                           Hugo_Symbol            = "gene",
                           Variant_Classification = "effect",
                           Variant_Type           = "Variant_Type")
#> Aggregation function missing: defaulting to length
#> Aggregation function missing: defaulting to length
#> Aggregation function missing: defaulting to length
#> Aggregation function missing: defaulting to length

# NOTE: IOBR provides mutation matrix(`tcga_stad_var`), if MAF data or UCSC can not be accessed. 
tcga_stad_var[1:5,1:5]
#>              ABCA12 ABCA13 ABCA2 ABCB1 ABCC9
#> TCGA-3M-AB46      0      0     0     0     0
#> TCGA-3M-AB47      0      0     0     0     0
#> TCGA-B7-5816      0      2     0     1     1
#> TCGA-B7-5818      0      0     0     0     0
#> TCGA-B7-A5TI      0      0     0     1     0


names(mut_list)
#> [1] "all"        "snp"        "indel"      "frameshift"
mut_list$all[1:5, 1:10]
#>              A1BG A1CF A2M A2ML1 A3GALT2 A4GALT A4GNT AAAS AACS AADAC
#> TCGA-3M-AB46    1    1   0     0       0      0     0    0    0     0
#> TCGA-3M-AB47    0    0   0     0       0      0     0    0    0     0
#> TCGA-B7-5816    0    0   1     0       0      0     0    1    0     0
#> TCGA-B7-5818    0    0   0     0       0      0     0    0    0     0
#> TCGA-B7-A5TI    0    0   0     0       0      0     0    0    0     0

# choose SNP mutation matrix as input
mut<-mut_list$snp
# NOTE: If the maximum of mutation counts of a single gene of a person in is over 4, the mutation data would be standardized according to following principles.
# mut[mut>=3&mut<=5]<-3
# mut[mut>5]<-4
# Each gene of every samples is categorized as binary variables(mutation or non-mutation) in the MAF data. 
##########################
res<-find_mutations(mutation_matrix     = mut, 
                    signature_matrix    = tcga_stad_sig,
                    id_signature_matrix = "ID",
                    signature           = "CD_8_T_effector",
                    min_mut_freq        = 0.01,
                    plot                = TRUE,
                    method              = "Wilcoxon",
                    save_path           = paste0("CD_8_T_effector-relevant-mutations"),
                    palette             = "jco",
                    show_plot           = T)
#> [1] ">>>> Result of Wilcoxon test"
#>             p.value  names statistic adjust_pvalue
#> PIK3CA 1.921035e-10 PIK3CA      4125  9.605174e-08
#> TCHH   1.961642e-05   TCHH      3312  4.904106e-03
#> SPEG   3.532750e-05   SPEG      1947  5.887916e-03
#> LRP1   7.511741e-05   LRP1      2649  9.389676e-03
#> WDFY3  1.257659e-04  WDFY3      2964  1.257659e-02
#> ARID1A 2.468609e-04 ARID1A      4878  2.057174e-02
#> PLXNA4 4.215809e-04 PLXNA4      3638  3.011292e-02
#> ANK3   6.399572e-04   ANK3      4446  3.933979e-02
#> DMD    7.364591e-04    DMD      5311  3.933979e-02
#> PLEC   8.026240e-04   PLEC      5562  3.933979e-02


if(T){
  library("IOBR")
  
  data("imvigor210_sig")
  data("imvigor210_pdata")
  # For analyses of binary variables
  input<-imvigor210_pdata %>% 
    dplyr::select(ID,BOR_binary) %>% 
    inner_join(.,imvigor210_sig,by="ID") %>% 
    filter(!is.na(.$BOR_binary)) %>% 
    filter(!.$BOR_binary=="NA")

  # Feature engineering
  res<-batch_wilcoxon(data    = as.data.frame(input),
                      target  = "BOR_binary",
                      group_names = c("NR","R"),
                      feature = colnames(input)[3:ncol(input)])
  head(res)
  #> # A tibble: 6 x 8
  #>   sig_names            p.value     NR      R statistic   p.adj log10pvalue stars
  #>   <chr>                  <dbl>  <dbl>  <dbl>     <dbl>   <dbl>       <dbl> <fct>
  #> 1 Mismatch_Repair      9.21e-6 -0.111  0.246    -0.357 0.00138        5.04 **** 
  #> 2 Cell_cycle           1.05e-5 -0.297  0.673    -0.969 0.00138        4.98 **** 
  #> 3 TMEscoreB_plus       1.35e-5  0.149 -0.432     0.581 0.00138        4.87 **** 
  #> 4 DNA_replication      1.38e-5 -0.147  0.353    -0.500 0.00138        4.86 **** 
  #> 5 DDR                  1.52e-5 -0.229  0.532    -0.760 0.00138        4.82 **** 
  #> 6 Nucleotide_excisio…  2.16e-5 -0.106  0.244    -0.350 0.00154        4.67 ****
  model_feas<-as.character(res[res$p.value<0.05,]$sig_names)
  
  input<-as.data.frame(imvigor210_sig)
  feas<-colnames(input)[colnames(input)%in%model_feas]
  input<-input[, c("ID",feas)]
  
  # target data
  pdata_group <- imvigor210_pdata[!imvigor210_pdata$BOR_binary=="NA",c("ID","BOR_binary")]
  pdata_group$BOR_binary <- ifelse(pdata_group$BOR_binary == "R", 1, 0)
  
  #Feature selection
  binomial_result <- BinomialModel(x           = input, 
                                   y           = pdata_group, 
                                   seed        = "123456", 
                                   scale       = TRUE,
                                   train_ratio = 0.7, 
                                   nfold       = 10, 
                                   plot        = T)
  
}


plot(binomial_result$lasso_result$model)

plot(binomial_result$ridge_result$model)


lapply(binomial_result[1:2], function(z)z$AUC)
#> $lasso_result
#>       lambda.min lambda.1se
#> train  0.8911641        0.5
#> test   0.7017974        0.5
#> 
#> $ridge_result
#>       lambda.min lambda.1se
#> train  0.7915115        0.5
#> test   0.7263072        0.5


data("imvigor210_pdata")
pdata_prog<-imvigor210_pdata %>% 
  dplyr::select(ID, OS_days, OS_status) %>% 
  filter(!is.na(.$OS_days)) %>% 
  filter(!.$OS_days=="NA") %>% 
  dplyr:: rename(time = OS_days) %>% 
  dplyr:: rename(status = OS_status) %>% 
  mutate(time = as.numeric(.$time)) %>%
  mutate(status = as.numeric(.$status))

pdata_prog<-as.data.frame(pdata_prog)
input<-as.data.frame(imvigor210_sig)

prognostic_result <- PrognosticModel(x           = input, 
                                     y           = pdata_prog, 
                                     scale       = T, 
                                     seed        = "123456", 
                                     train_ratio = 0.8,
                                     nfold       = 10,
                                     plot        = T)

plot(prognostic_result$lasso_result$model)


#######https://mp.weixin.qq.com/s?__biz=MzUyODk4MTIyNA==&mid=2247484124&idx=1&sn=e827fd4f71fab3cae14d7aae70ddacec&chksm=fa694f5bcd1ec64d82d6a875b9a2469d83ed3f253f27fae2b0dd19bbcaebff7f8d6630536e67&scene=178&cur_album_id=1661319977951346691#rd
if(F){
library(UCSCXenaTools)
eset_stad<-XenaGenerate(subset = XenaCohorts =="GDC TCGA Stomach Cancer (STAD)") %>% 
  XenaFilter(filterDatasets    = "TCGA-STAD.htseq_counts.tsv") %>% 
  XenaQuery() %>%
  XenaDownload() %>% 
  XenaPrepare()

# Remove the version numbers in Ensembl ID.
eset_stad$Ensembl_ID<-substring(eset_stad$Ensembl_ID, 1, 15)
eset_stad<-column_to_rownames(eset_stad, var = "Ensembl_ID")
# Revert back to original format because the data from UCSC was log2(x+1)transformed.
eset_stad<-(2^eset_stad)+1

eset_stad<-count2tpm(countMat = eset_stad, idType = "Ensembl", source = "default")
# TCGA-D7-5577-01A TCGA-D7-6818-01A TCGA-BR-7958-01A TCGA-D7-8572-01A TCGA-VQ-A91Z-01A
# RPS20         72997.5064        42582.747        53479.686        50017.591       56903.9292
# SLC25A5       50356.7131        53582.503        42953.670        40360.665       47658.4941
# HSPB6           653.1407         5542.878         4310.936         8068.731         502.1561
# GPRC5A        20351.1750         1668.879        16001.004        20881.297        9891.7909
# CSDE1          8552.9151        11528.998        14264.571        20035.149       14262.3733


library("IOBR")
signature_metabolism[1:3]

if(T){
# 1. 所以只需要构建list即可，如果你的signature gene set只有两三个，可以直接将他们构建一个list对象:
my_signature<-list("CD8" = c( "CD8A", "GZMA", "GZMB", "IFNG", "CXCL9", "CXCL10", "PRF1", "TBX21"),
                   "Immune_Checkpoint" = c("CD274", "PDCD1LG2", "CTLA4", "PDCD1", "LAG3", "HAVCR2", "TIGIT"))
my_signature
}
if(T){
# 2. 或者直接将signature添加到IOBR已有的signature中一起计算：

#别随便修改谢谢
# signature_collection$CD8<-c( "CD8A", "GZMA", "GZMB", "IFNG", "CXCL9", "CXCL10", "PRF1", "TBX21")
}

data(eset_stad)
}
# my_signature=signature_collection[1:5]
score<-calculate_sig_score(eset            = genes_expr, 
                           signature       = my_signature,
                           method          = "integration",
                           mini_gene_count = 2)
?calculate_sig_score
score[1:5,1:5]
# A tibble: 5 x 5
#   Index ID           Nature_metabolism_Hypoxia Winter_hypoxia_signature Hu_hypoxia_signature
# <int> <fct>                            <dbl>                    <dbl>                <dbl>
  # 1     1 TCGA-B7-5818                     1.65                    1.25                -1.16  
  # 2     2 TCGA-BR-4187                    -0.218                   0.0715               0.0982
  # 3     3 TCGA-BR-4201                     1.18                    0.294                1.28  
  # 4     4 TCGA-BR-4253                     0.642                  -1.16                -1.53  
  # 5     5 TCGA-BR-4256                     0.528     
  # 加载R包
  if (!requireNamespace("IOBR", quietly = TRUE)) 
    devtools::install_github("IOBR/IOBR")
library(IOBR)
library(dplyr)

#构建文件夹用于储存解析文件
#####################################
path<-"1-TME-signature-deconvolution"
if (!file.exists(path)) dir.create(path)
abspath<-paste0(getwd(),"/" ,path, "/")

#读入标准化好的基因表达矩阵: 行名为人的基因symbol，列为样本
# 如何获取矩阵请参考既往链接
eset<-as.matrix(genes_expr)
eset[1:5, 1:5]
#           TCGA-B7-5818 TCGA-BR-4187 TCGA-BR-4201 TCGA-BR-4253 TCGA-BR-4256
# MT-CO1      15.18012     15.55806     14.60960     14.63728     15.23528
# MT-CO3      14.75536     15.19199     14.55337     13.54925     14.30425
# MT-ND4      14.19637     15.61564     15.80262     14.98329     14.83764
# MT-CO2      15.10790     15.50514     15.43261     14.52009     14.65806
# MT-RNR2     14.22690     14.71157     13.48096     13.32553     13.55689

####################################
array<-FALSE  #是否是芯片的数据，如果是RNAseq, 此处为FALSE
ProjectID<-"TCGA-COAD"
tumor_type<-"coad" #命名癌种用于TIMER的解析

#各种微环境工具的解析
#######################################
cibersort<-deconvo_tme(eset = eset,method = "cibersort",arrays = array,perm = 1000 )
epic     <-deconvo_tme(eset = eset,method = "epic",arrays = array)
mcp      <-deconvo_tme(eset = eset,method = "mcpcounter")
xcell    <-deconvo_tme(eset = eset,method = "xcell",arrays = array)
estimate <-deconvo_tme(eset = eset,method = "estimate")
timer    <-deconvo_tme(eset = eset,method = "timer",group_list = rep(tumor_type,dim(eset)[2]))
quantiseq<-deconvo_tme(eset = eset,method = "quantiseq", tumor = TRUE, arrays = array, scale_mrna = TRUE)
ips      <-deconvo_tme(eset = eset,method = "ips",plot= FALSE)

#合并数据
tme_combine<-cibersort %>% 
  inner_join(.,mcp,by       = "ID") %>% 
  inner_join(.,xcell,by     = "ID") %>%
  inner_join(.,epic,by      = "ID") %>% 
  inner_join(.,estimate,by  = "ID") %>% 
  inner_join(.,quantiseq,by = "ID") %>% 
  inner_join(.,timer,by     = "ID") %>% 
  inner_join(.,ips,by       = "ID") 

# tme_combine<-tme_combine[,-c(grep(colnames(tme_combine),pattern = "Index"))]
############################################
save(tme_combine,file = paste0(abspath,"1-",ProjectID,"-TME-Cell-fration.RData"))
print(paste0( ">>>>> TME cell deconvolution was finished: ", ProjectID))


#使用多个方法计算signature-score
#######################################
help("calculate_sig_score")
sig_res<-calculate_sig_score(pdata = NULL,
                             eset = eset,
                             signature = signature_collection,
                             method = "integration",
                             adjust_eset = T,
                             mini_gene_count = 2)
print(paste0( ">>>>> Signatures score esitmation of IOBR collection was finished: ", ProjectID))
save(sig_res,file = paste0(abspath,"2-",ProjectID,"-Signature-score-mycollection.RData"))
#使用ssgsea的方法计算MsigDb的signature-score
########################################
sig_go_kegg<-calculate_sig_score(pdata = NULL,
                                 eset = eset,
                                 signature = c(hallmark,go_bp,go_cc,go_mf,kegg,reactome),
                                 method = "ssgsea",
                                 mini_gene_count = 2)


save(sig_go_kegg,file = paste0(abspath,"3-",ProjectID,"-Signature-score-Hallmark-GO-KEGG.RData"))
print(paste0( ">>>>> HALLMARK GO KEGG REACTOME esitmation were finished: ", ProjectID))
#########################################

#合并所有的解析数据
tme_sig_combin<-tme_combine %>% 
  inner_join(.,sig_res,by = "ID") %>% 
  inner_join(.,sig_go_kegg,by = "ID")
save(tme_sig_combin,file = paste0(abspath,"0-",ProjectID,"-Merge-TME-Signature-and-Hallmark-GO-KEGG.RData"))