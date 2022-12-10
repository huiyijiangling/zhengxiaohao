# Step2 clinical information --------------------------

Rdata_file <- paste('./data/', tumor_type, '.phenoData.Rdata', sep = '')
if (!file.exists( Rdata_file )) {
  phenoData <- read.table( destfile,
                           header = T,
                           sep = '\t',
                           quote = '' )
  rownames( phenoData ) <- phenoData[ , 1]
  ## name <- gsub(pattern = "-", replacement = ".", name)
  colnames( phenoData )[1] <- "Tumor_Sample_Barcode"
  phenoData[1:5, 1:5]
  save( phenoData, file = Rdata_file )
}else{
  load( Rdata_file )
}


回归分析前数据准备


# step1 Keep only tumor samples -------------------------------------------

ncol(raw_data)
group_list <- factor(
  ifelse( as.numeric( substr( 
    colnames( raw_data ),14,15)) < 10, 'tumor', 'normal'))
table(group_list)
AssayData <- na.omit(raw_data)

AssayData <- AssayData[, group_list == 'tumor']
dim(AssayData)
AssayData[1:3, 1:3]



# step2 Only retain lncRNA ------------------------------------------------

## ENSEMBLTO SYMBOL
library( "clusterProfiler" )
library( "org.Hs.eg.db" )
keytypes(org.Hs.eg.db)
library("stringr")
rownames( AssayData ) <- str_sub(rownames( AssayData ), start = 1, end = 15)
AssayData$ENSEMBL <- rownames( AssayData )
gene_name <- bitr( AssayData$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", 
                   OrgDb = org.Hs.eg.db )
head( gene_name )
AssayData <- AssayData[gene_name$ENSEMBL, ]
dim(AssayData)
AssayData[1:3, 1:3]

## Pick the gene with the highest expression
gene_name$max <- apply(AssayData, 1, max)
gene_name <- gene_name[order(gene_name$SYMBOL,
                             gene_name$max,
                             decreasing = T), ]
gene_name <- gene_name[!duplicated(gene_name$SYMBOL), ]
dim( gene_name )

## pick lncRNA
{
  gene2type = read.table( './raw_data/gencode.v25lift37.annotation.gtf.gene2type' )
  colnames( gene2type ) = c( "gene", "type" )
}

dim( gene2type )
sort( table( gene2type$type ) )
gene2type = gene2type[ gene2type[,2] == 'lincRNA', ]
length( unique( gene2type$gene ) )

gene_name <- gene_name[gene_name$SYMBOL %in% gene2type$gene, ]
length(gene_name$SYMBOL)

AssayData <- AssayData[gene_name$ENSEMBL, ]

rownames(AssayData) <- gene_name$SYMBOL



# step3 Random grouping ---------------------------------------------------

set.seed(566)
k <- sample(1:ncol(AssayData), ncol(AssayData) / 2 )

t_exp <- AssayData[,k]
v_exp <- AssayData[,-k]



## Based on a data set: t_exp, v_exp, AssayData


# step4 Extract clinical information related to risk factors --------------

pheno <- phenoData[colnames(t_exp), ]
head(pheno)
colnames(pheno)
pheno[1:10, 137:138]

## race
pheno$race <- str_split(pheno$race, ' ', simplify = T)[, 1]
table(pheno$race)

## OS.time
pheno[, 137][is.na(pheno[,137])] = 0
pheno[, 138][is.na(pheno[, 138])] = 0
pheno$OS.days <- as.numeric(pheno[, 137]) + as.numeric(pheno[, 138])
pheno$time <- pheno$OS.days / 30
boxplot(pheno$time)

## status
pheno$status <- ifelse(pheno$vital_status.diagnoses == "alive", 0, 1)

## age
pheno$year_of_death.demographic[is.na(pheno$year_of_death.demographic)] <- '2016'
table(pheno$year_of_death.demographic)
pheno$age <- as.numeric(pheno$year_of_death.demographic) - as.numeric(pheno$year_of_birth.demographic)
boxplot(pheno$age)
median_age <- sort(pheno$age)[length(pheno$age) / 2]
pheno$age_group <- ifelse( pheno$age >= median_age, 'older', 'younger' )
table(pheno$age_group)

## stage
library(stringr) 
pheno$tumor_stage <- str_split(pheno$tumor_stage.diagnoses, ' ', simplify = T)[, 2]
pheno$tumor_stage <- ifelse(pheno$tumor_stage %in% c("i", "ii"),
                            "early_stage", "later_stage")
table(pheno$tumor_stage)

## bind your clinical data
pheno <- pheno[,c("Tumor_Sample_Barcode",
                  "race",
                  "gender.demographic",
                  "tumor_stage",
                  "status",
                  "time",
                  "age_group")]
colnames(pheno) <- c('ID', 'race', 'gender', 'tumor_stage', 'status', 'time', 'age_group')

pheno$ID <- toupper(pheno$ID)
dim(pheno)
head(pheno)
t_exp[1:4, 1:4]


COX回归分析


首先挑选出感兴趣的基因



# step5 COX analysis ------------------------------------------------------
## http://www.sthda.com/english/wiki/cox-proportional-hazards-model

library(survival)
## Eliminate more than half of the gene data not expressed
t_exp <- t_exp[apply(t_exp, 1, function(x) sum(x == 0)) < (ncol(t_exp)*0.49), ]
dim(t_exp)
## [1] 445 187

## single
group_data <- apply(t_exp , 1 , function(gene){
  name <- rownames(gene)
  gene <- unlist(gene)
  group <- ifelse(gene >= median(gene), 'high', 'low')
  names(group) <- name
  return(group)
})
group_data <- as.data.frame(group_data, stringsAsFactors = F)
survival_dat <- data.frame(race = pheno$race,
                           gender = pheno$gender,
                           stage = pheno$tumor_stage, 
                           age_group = pheno$age_group,
                           status = pheno$status,
                           time = pheno$time,
                           stringsAsFactors = F)
survival_dat <- cbind(group_data, survival_dat)
colnames(survival_dat) <- sub("\\-", "", colnames(survival_dat))
covariates <- as.character(colnames(survival_dat))
univ_formulas <- sapply(covariates,
                        function(x){
                          ##print(x)
                          as.formula(paste('Surv(time, status)~', x))
                        })

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = survival_dat)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value <- signif(x$wald["pvalue"], digits = 2)
                         beta <- signif(x$coef[1], digits = 2)
                         HR <- signif(x$coef[2], digits = 2)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res <- c(beta, HR, p.value)
                         names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
                         return(res)
                       })
res_single <- as.data.frame(t(do.call(cbind, univ_results)))
table(res_single$p.value <= 0.01)
res_single <- res_single[res_single$p.value <= 0.01, ]
res_single <- res_single[order(res_single$p.value), ]
single_pick <- rownames(res_single)

## multi
multi_results <- apply(t_exp , 1 , function(gene){
  ## gene <- t_exp[1, ]
  gene <- unlist(gene)
  group <- ifelse(gene >= median(gene), 'high', 'low')
  survival_dat <- data.frame(group = group,
                             race = pheno$race,
                             gender = pheno$gender,
                             stage = pheno$tumor_stage, 
                             age_group = pheno$age_group,
                             status = pheno$status,
                             time = pheno$time,
                             stringsAsFactors = F)
  res.cox <- coxph(Surv(time, status) ~ age_group + stage + gender + race + group, 
                   data =  survival_dat)
  ## summary(res.cox)
  beta <- coef(res.cox)
  se <- sqrt(diag(vcov(res.cox)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  res <- as.data.frame(round(cbind(coef = beta,
                                   se = se,
                                   z = beta/se,
                                   p.value = 1 - pchisq((beta/se)^2, 1),
                                   HR = HR,
                                   HRse = HRse,
                                   HRz = (HR - 1) / HRse,
                                   HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                                   HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                                   HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 5))
  return(res['grouplow',])
})
multi_results <- do.call(rbind, multi_results)
table(multi_results$p.value <= 0.01)
res_multi <- multi_results[multi_results$p.value <= 0.01, ]
res_multi <- res_multi[order(res_multi$p.value), ]
multi_pick <- rownames(res_multi)

overgene <- intersect(multi_pick, single_pick)


森林图
model_exp <- t(log2(t_exp[overgene,] + 1))
colnames(model_exp) <- overgene
dat <- cbind(pheno, model_exp)

dat$gender <- factor(dat$gender)
dat$tumor_stage <- factor(dat$tumor_stage)

library("survminer")
colnames(dat)
s <- Surv(time, status) ~ race + gender + tumor_stage + age_group + LINC01136 + LINC00346 + C2orf48 + LINC01139 + LINC01559
s <- Surv(time, status) ~ LINC01136 + LINC00346 + C2orf48 + LINC01139 + LINC01559

model <- coxph(s, data = dat )
summary(model, data = dat)
options(scipen = 1)
ggforest(model, data = dat, 
         main = "Hazard ratio", 
         cpositions = c(0.10, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)


生存曲线图
model_exp <- t_exp[overgene, ]
for (i in 1:nrow(model_exp)) {
  gene <- model_exp[i, ]
  name <- rownames(gene)
  gene <- unlist(gene)
  group <- ifelse(gene >= median(gene), 'high', 'low')
  survival_dat <- data.frame(group = group,
                             status = pheno$status,
                             time = pheno$time,
                             stringsAsFactors = F)
  fit <- survfit(Surv(time, status) ~ group, data = survival_dat)
  ggsurvplot(fit, data = survival_dat,
             surv.median.line = "hv",
             legend.title = "Group",
             legend.labs = c("High", "Low"),
             pval = TRUE,
             conf.int = TRUE,
             palette = "jco",
             ggtheme = theme_bw()
  )
  ggsave(filename = paste('./fig/', tumor_type, name, '.png', sep = ''))
}




ROC曲线
library(lars) 
library(glmnet) 
x <- t(log2(t_exp[overgene,] + 1))
y <- survival_dat$status 
cv_fit <- cv.glmnet(x = x, y = y, alpha = 1, nlambda = 1000) 
lasso.prob <- predict(cv_fit,
                      newx = x , 
                      s = c(cv_fit$lambda.min, cv_fit$lambda.1se) )


library(timeROC)
pheno$fp <- as.numeric(lasso.prob[,1])
with(pheno,
     ROC <<- timeROC(T = time,
                     delta = status,
                     marker = fp,
                     cause = 1,
                     weighting = "marginal",
                     times = c(20, 60),
                     ROC = TRUE,
                     iid = TRUE)
)

plot(ROC, time = 60, col = "blue", add = F)
plot(ROC, time = 20, col = "red", add = T)