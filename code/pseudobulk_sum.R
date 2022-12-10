# all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
#                                                      n_genes = n_genes,
#                                                      n_per_group = n_per_group,
#                                                      n_cases = n_cases,
#                                                      n_controls = n_controls,
#                                                      cells_per_case = cells_per_case,
#                                                      cells_per_control = cells_per_control,
#                                                      ncells_variation_type = ncells_variation_type,
#                                                      foldchange = foldchange,
#                                                      decrease_dropout = decrease_dropout,
#                                                      alter_dropout_cases = alter_dropout_cases))
# 
# genecounts <- as.matrix(all_genes[,c(-1,-2,-3)])
# genecounts <- genecounts[ ,which(apply(genecounts, 2, mean) > 5)]
# genecounts <- cbind(all_genes[,1:2],genecounts)
# 
# computesums <- function(x){tapply(x,genecounts[,2],sum)}
# cellsums <- sapply(genecounts[,c(-1,-2)],computesums)
# rownames(cellsums) <- paste0(rownames(cellsums),"_sum")
# coldata <- as.data.frame(cbind(rownames(cellsums),rownames(cellsums)))
# colnames(coldata) <- c("SampleID","ToSep")
# coldata <- tidyr::separate(coldata,ToSep,c("Status", "Donor_Number", "sum"), sep="_")
# rownames(coldata) <- coldata$SampleID
# coldata$Status <- as.factor(coldata$Status)
# coldata$Status <- stats::relevel(coldata$Status, "Control")
# cellsums <- round(t(cellsums),0)
# cellsums <- cellsums[, rownames(coldata)]
# 
# #
# all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
#                                                      n_genes = n_genes,
#                                                      n_per_group = n_per_group,
#                                                      n_cases = n_cases,
#                                                      n_controls = n_controls,
#                                                      cells_per_case = cells_per_case,
#                                                      cells_per_control = cells_per_control,
#                                                      ncells_variation_type = ncells_variation_type,
#                                                      foldchange = foldchange,
#                                                      decrease_dropout = decrease_dropout,
#                                                      alter_dropout_cases = alter_dropout_cases))
# 
# genecounts <- as.matrix(all_genes[,c(-1,-2,-3)])
# genecounts <- genecounts[ ,which(apply(genecounts, 2, mean) > 5)]
# genecounts <- cbind(all_genes[,1:2],genecounts)
# 
# computemeans <- function(x){tapply(x,genecounts[,2],mean)}
# cellmeans <- sapply(genecounts[,c(-1,-2)],computemeans)
# rownames(cellmeans) <- paste0(rownames(cellmeans),"_Mean")
# coldata <- as.data.frame(cbind(rownames(cellmeans),rownames(cellmeans)))
# colnames(coldata) <- c("SampleID","ToSep")
# coldata <- tidyr::separate(coldata,ToSep,c("Status", "Donor_Number", "Mean"), sep="_")
# rownames(coldata) <- coldata$SampleID
# coldata$Status <- as.factor(coldata$Status)
# coldata$Status <- stats::relevel(coldata$Status, "Control")
# cellmeans <- round(t(cellmeans),0)
# cellmeans <- cellmeans[, rownames(coldata)]





# # README
# 
# Libra is an R package to perform differential expression on single-cell data. Libra implements a total of 22 unique differential expression methods that can all be accessed from one function. These methods encompass traditional single-cell methods as well as methods accounting for biological replicate including pseudobulk and mixed model methods. The code for this package has been largely inspired by the [Seurat](https://satijalab.org/seurat/) and [Muscat](https://github.com/HelenaLC/muscat) packages. Please see the documentation of these packages for further information.
# 
# 
# 
# ## Usage
# 
# The main function of Libra, `run_de`, takes as input a preprocessed features-by-cells (e.g., genes-by-cells for scRNA-seq) matrix, and a data frame containing metadata associated with each cell, minimally including the cell type annotations, replicates, and sample labels to be predicted.
# This means that in order to use Libra, you should have pre-processed your data (e.g., by read alignment and cell type assignment for scRNA-seq) across all experimental conditions.
# 
# Libra provides a universal interface to perform differential expression using 22 discrete methods. These methods are summarized as follows:
#   
#   __Single cell methods__
# 
# - Wilcoxon Rank-Sum test
# - Likelihood ratio test
# - Student's t-test
# - Negative binomial linear model
# - Logistic regression
# - MAST
# 
# __Pseudobulk methods__
# 
# - edgeR-LRT
# - edgeR-QLF
# - DESeq2-LRT
# - DESeq2-Wald
# - limma-trend
# - limma-voom
# 
# __Mixed model methods__
# 
# - Linear mixed model
# - Linear mixed model-LRT
# - Negative binomial generalized linear mixed model
# - Negative binomial generalized linear mixed model-LRT
# - Negative binomial generalized linear mixed model with offset
# - Negative binomial generalized linear mixed model with offset-LRT
# - Poisson generalized linear mixed model
# - Poisson generalized linear mixed model-LRT
# - Poisson generalized linear mixed model with offset
# - Poisson generalized linear mixed model with offset-LRT
# 
# By default Libra will use a pseudobulk approach, implementing the `edgeR` package with a likelihood ratio test (LRT) null hypothesis testing framework. Each of the 22 tests can be accessed through three key variables of the `run_de` function: `de_family`, `de_method`, and `de_type`. Their precise access arguments are summarized in the below table.
# 
# | Method                                                       | de_family  | de_method       | de_type |
# | ------------------------------------------------------------ | ---------- | --------------- | ------- |
# | Wilcoxon Rank-Sum test                                       | singlecell | wilcox          |         |
# | Likelihood ratio test                                        | singlecell | bimod           |         |
# | Student's t-test                                             | singlecell | t               |         |
# | Negative binomial linear model                               | singlecell | negbinom        |         |
# | Logistic regression                                          | singlecell | LR              |         |
# | MAST                                                         | singlecell | MAST            |         |
# | edgeR-LRT                                                    | pseudobulk | edgeR           | LRT     |
# | edgeR-QLF                                                    | pseudobulk | edgeR           | QLF     |
# | DESeq2-LRT                                                   | pseudobulk | DESeq2          | LRT     |
# | DESeq2-Wald                                                  | pseudobulk | DESeq2          | Wald    |
# | limma-trend                                                  | pseudobulk | limma           | trend   |
# | limma-voom                                                   | pseudobulk | limma           | voom    |
# | Linear mixed model                                           | mixedmodel | linear          | Wald    |
# | Linear mixed model-LRT                                       | mixedmodel | linear          | LRT     |
# | Negative binomial generalized linear mixed model             | mixedmodel | negbinom        | Wald    |
# | Negative binomial generalized linear mixed model-LRT         | mixedmodel | negbinom        | LRT     |
# | Negative binomial generalized linear mixed model with offset | mixedmodel | negbinom_offset | Wald    |
# | Negative binomial generalized linear mixed model with offset-LRT | mixedmodel | negbinom_offset | LRT     |
# | Poisson generalized linear mixed model                       | mixedmodel | poisson         | Wald    |
# | Poisson generalized linear mixed model-LRT                   | mixedmodel | poisson         | LRT     |
# | Poisson generalized linear mixed model with offset           | mixedmodel | poisson_offset  | Wald    |
# | Poisson generalized linear mixed model with offset-LRT       | mixedmodel | poisson_offset  | LRT     |
#   
#   If batch effects are present in the data, these should be accounted for, e.g., using [Seurat](https://www.sciencedirect.com/science/article/pii/S0092867419305598) or [Harmony](https://www.nature.com/articles/s41592-019-0619-0), to avoid biasing differential expression by technical differences or batch effects.
# 
# To run Libra with default parameters on a genes-by-cells scRNA-seq matrix `expr`, and an accompanying data frame `meta`, with `cell_type`, `replicate`, and `label` columns containing cell types, replicates, and experimental conditions, respectively, use the `run_de` function:
#   
# replicate_col = "replicate",
# cell_type_col = "cell_type",
# label_col = "label",

#to_pseudobulk 用的是sum 
#Al-Murphy/reanalysis_scRNA_seq_benchmark 说 to_pseudobulk mean 更好 
# 但我还是觉得sum 更符合原理
library(Libra)
#pseudobulk | edgeR           | LRT 
#pseudobulk | limma           | voom
#singlecell | wilcox 

data("hagai_toy")

head(hagai_toy@meta.data)

# 如果您的列有不同的名称，您可以使用 、 和 参数指定cell_type_col这些replicate_col名称label_col：
# 
# > DE = run_de(expr, meta = meta, cell_type_col = "cell.type", label_col = "condition")
# 如果您想将伪散装矩阵存储在变量中，在运行微分表达式之前，您可以执行以下操作：
# 
# > matrices = to_pseudobulk(expr, meta = meta)

DE = run_de(hagai_toy)
head(DE)
#
DE = run_de(hagai_toy, de_family = 'pseudobulk', de_method = 'DESeq2', de_type = 'LRT', n_threads = 16)
head(DE)
```

Alternatively, we can use a mixed-model approach, which by default will use a negative binomial model structure:
  
  ```r
> DE = run_de(hagai_toy, de_family = 'mixedmodel')
> head(DE)
```

However, this can be adjusted using the `de_method` argument:
  
  ```r
> DE = run_de(hagai_toy, de_family = 'mixedmodel', de_method = 'linear', de_type = 'LRT', n_threads = 16)
> head(DE)
```

Running this example on a MacBook should be instantaneous.
However, analyzing >20 real single-cell RNA-seq datasets, we found Libra takes a median of ~5 minutes.
In general, runtime scales close to linearly with the number of cell _types_ and _cells_.
If using _mixed models_, by default, Libra runs on four cores, with each gene analyzed on a different core.
To change the number of cores, use the `n_threads` argument.
For example, running Libra on eight threads:
  

DE = run_de(hagai_toy, de_family = 'mixedmodel', de_method = 'linear', de_type = 'LRT', n_threads = 8)

计算增量方差
我们最近表明，差异表达的统计方法必须考虑生物复制的内在变异性，才能在单细胞数据中产生生物学上准确的结果（Squair 等人，2021，Biorxiv；https: //www.biorxiv.org/content/ 10.1101/2021.03.12.435024v1）。在相同的实验条件下，复制表现出基因表达的内在差异，这反映了生物学和技术因素。我们推断，未能解释这些差异可能会导致方法错误地将复制之间的固有变异性归因于扰动的影响。为了研究这种可能性，我们比较了假体和假复制中每个基因的表达差异。我们将此度量称为“增量方差”。用户可以使用 delta 方差的计算来告知他们的差异表达结果。例如，通过不考虑生物复制的方法（即“单细胞”方法）鉴定为差异表达的基因具有高 delta 方差，应谨慎处理，因为它们可能是假阳性。

 DV = calculate_delta_variance(hagai_toy)
该函数将返回一个向量列表，每个细胞类型对应一个向量，每个向量都包含输入表达矩阵中存在的基因的 delta 方差。

元

随附的元数据，其中行名称与输入的列名匹配。

to_pseudobulk(
  input,
  meta = NULL,
  replicate_col = "replicate",#生物学复制
  cell_type_col = "cell_type",#细胞类型
  label_col = "label",#实验标签
  min_cells = 3,
  min_reps = 2,
  min_features = 0#保留基因的最小计数数。默认为0
)
lapply(hagai_toy@meta.data, table)

# 
# $nCount_RNA
# 
# 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 
# 2  7  5 14 18 23 29 29 27 24 35 31 27 24 38 25 30 20 31 23 22 15 21 17 16 11  4 12  4  4  2  3  1  2  1  1  1  1 
# 
# $nFeature_RNA
# 
# 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 19 
# 2  7  8 28 36 52 50 67 45 68 56 47 55 35 24  9  8  2  1 
# 
# $replicate
# 
# mouse1 mouse2 mouse3 
# 200    200    200 
# 
# $label
# 
# lps4 unst 
# 300  300 
# 
# $cell_type
# 
# bone marrow derived mononuclear phagocytes 
# 600 

run_de(
  input,
  meta = NULL,
  replicate_col = "replicate",
  cell_type_col = "cell_type",
  label_col = "label",
  min_cells = 3,
  min_reps = 2,
  min_features = 0,
  de_family = "pseudobulk",
  de_method = "edgeR",
  de_type = "LRT",
  n_threads = 2
)
