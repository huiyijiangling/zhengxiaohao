devtools::install_github("Moonerss/CIBERSORT")
#直接编程
library(CIBERSORT)
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
results <- cibersort(sig_matrix, mixture_file)
## example 2

data(LM22)
data(mixed_expr)
results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr,perm=0,QN=T)
# cibersort 是可以内部 QN的