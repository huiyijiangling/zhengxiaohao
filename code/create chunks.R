n = 3
split(x, sort(x%%n))
library(parallel)
splitIndices(20, 3)
split(x, rep_len(1:n, length(x)))
split(x, sort(rep_len(1:n, length(x))))