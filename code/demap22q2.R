#depmap 写入
sample_info_21Q4 <- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\21Q4\sample_info.csv)",header = T,sep = ",",verbose = F, integer64 = 'numeric')
sample_info_21Q4 = as.data.frame(sample_info_21Q4)
sample_info_21Q4[sample_info_21Q4=="Unknown"] <- NA
sample_info_21Q4[sample_info_21Q4=="unknown"] <- NA
sample_info_21Q4[sample_info_21Q4==""] <- NA

C:\Users\zxh\Desktop\R\depmap\22Q2\DOWNLOAD FROM FIRST PAGE\CRISPR_(DepMap_22Q2_Public,_Chronos).csv


sample_info_21Q4 <- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\DOWNLOAD FROM FIRST PAGE\CRISPR_(DepMap_22Q2_Public,_Chronos).csv)",header = T,sep = ",",verbose = F, integer64 = 'numeric')
sample_info_21Q42 <- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\DOWNLOAD FROM FIRST PAGE\CRISPR_(DepMap_22Q2_Public+Score,_Chronos).csv)",header = T,sep = ",",verbose = F, integer64 = 'numeric')
sample_info_21Q4 <- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\CRISPR_gene_effect.csv)",header = T,sep = ",",verbose = F, integer64 = 'numeric')
sample_info_21Q42 <- data.table::fread(file =r"(C:\Users\zxh\Desktop\R\depmap\22Q2\CRISPR_gene_dependency.csv)",header = T,sep = ",",verbose = F, integer64 = 'numeric')


sample_info_21Q4 = as.data.frame(sample_info_21Q4)
sample_info_21Q4[sample_info_21Q4=="Unknown"] <- NA
sample_info_21Q4[sample_info_21Q4=="unknown"] <- NA
sample_info_21Q4[sample_info_21Q4==""] <- NA
#957
#1086
datasetid=read.table("h:/gsedownload.txt",header=T)
