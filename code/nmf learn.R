library(NMF)
# generate a synthetic dataset with known classes: 50 features, 23 samples (10+5+8)
n <- 20; counts <- c(5, 3, 2);
p <- sum(counts)
x <- syntheticNMF(n, counts)
dim(x)

# build the true cluster membership
groups <- unlist(mapply(rep, seq(counts), counts))

# run on a data.frame
res <- nmf(data.frame(x), 3)

# missing method: use algorithm suitable for seed
res <- nmf(x, 2, seed=rnmf(2, x))
algorithm(res)
res <- nmf(x, 2, seed=rnmf(2, x, model='NMFns'))
algorithm(res)

# compare some NMF algorithms (tracking the approximation error)
res <- nmf(x, 2, list('brunet', 'lee', 'nsNMF'), .options='t')
res
summary(res, class=groups)

# plot the track of the residual errors
plot(res)

# specify algorithm by its name
res <- nmf(x, 3, 'nsNMF', seed=123) # nonsmooth NMF
# names are partially matched so this also works
identical(res, nmf(x, 3, 'ns', seed=123))

res <- nmf(x, 3, 'offset') # NMF with offset

# run a custom algorithm defined as a standard function
myfun <- function(x, start, alpha){
  # update starting point
  # ...
  basis(start) <- 3 * basis(start)
  # return updated point
  start
}

res <- nmf(x, 2, myfun, alpha=3)
algorithm(res)
# error: alpha missing
try( nmf(x, 2, myfun) )

# possibly the algorithm fits a non-standard NMF model, e.g. NMFns model
res <- nmf(x, 2, myfun, alpha=3, model='NMFns')
modelname(res)

# assume a known NMF model compatible with the matrix `x`
y <- rnmf(3, x)
# fits an NMF model (with default method) on some data using y as a starting point
res <- nmf(x, y)
# the fit can be reproduced using the same starting point
nmf.equal(nmf(x, y), res)

# missing method: use default algorithm
res <- nmf(x, 3)

# Fit a 3-rank model providing an initial value for the basis matrix
nmf(x, rmatrix(nrow(x), 3), 'snmf/r')

# Fit a 3-rank model providing an initial value for the mixture coefficient matrix
nmf(x, rmatrix(3, ncol(x)), 'snmf/l')

# default fit
res <- nmf(x, 2)
summary(res, class=groups)

# run default algorithm multiple times (only keep the best fit)
res <- nmf(x, 3, nrun=10)
res
summary(res, class=groups)

# run default algorithm multiple times keeping all the fits
res <- nmf(x, 3, nrun=10, .options='k')
res
summary(res, class=groups)

## Note: one could have equivalently done
# res <- nmf(V, 3, nrun=10, .options=list(keep.all=TRUE))

# use a method that fit different model
res <- nmf(x, 2, 'nsNMF')
fit(res)

# pass parameter theta to the model via `...`
res <- nmf(x, 2, 'nsNMF', theta=0.2)
fit(res)

## handling arguments in `...` and model parameters
myfun <- function(x, start, theta=100){ cat("theta in myfun=", theta, "\n\n"); start }
# no conflict: default theta
fit( nmf(x, 2, myfun) )
# no conlfict: theta is passed to the algorithm
fit( nmf(x, 2, myfun, theta=1) )
# conflict: theta is used as model parameter
fit( nmf(x, 2, myfun, model='NMFns', theta=0.1) )
# conflict solved: can pass different theta to model and algorithm
fit( nmf(x, 2, myfun, model=list('NMFns', theta=0.1), theta=5) )

## USING SEEDING METHODS

# run default algorithm with the Non-negative Double SVD seeding method ('nndsvd')
res <- nmf(x, 3, seed='nndsvd')

## Note: partial match also works
identical(res, nmf(x, 3, seed='nn'))

# run nsNMF algorithm, fixing the seed of the random number generator
res <- nmf(x, 3, 'nsNMF', seed=123456)
nmf.equal(nmf(x, 3, 'nsNMF', seed=123456), res)

# run default algorithm specifying the starting point following the NMF standard model
start.std <- nmfModel(W=matrix(0.5, n, 3), H=matrix(0.2, 3, p))
nmf(x, start.std)

# to run nsNMF algorithm with an explicit starting point, this one
# needs to follow the 'NMFns' model:
start.ns <- nmfModel(model='NMFns', W=matrix(0.5, n, 3), H=matrix(0.2, 3, p))
nmf(x, start.ns)
# Note: the method name does not need to be specified as it is infered from the
# when there is only one algorithm defined for the model.

# if the model is not appropriate (as defined by the algorihtm) an error is thrown
# [cf. the standard model doesn't include a smoothing parameter used in nsNMF]
try( nmf(x, start.std, method='nsNMF') )

## Callback functions
# Pass a callback function to only save summary measure of each run
res <- nmf(x, 3, nrun=3, .callback=summary)
# the callback results are simplified into a matrix
res$.callback
res <- nmf(x, 3, nrun=3, .callback=summary, .opt='-S')
# the callback results are simplified into a matrix
res$.callback

# Pass a custom callback function
cb <- function(obj, i){ if( i %% 2 ) sparseness(obj) >= 0.5 }
res <- nmf(x, 3, nrun=3, .callback=cb)
res$.callback

# Passs a callback function which throws an error
cb <- function(){ i<-0; function(object){ i <<- i+1; if( i == 1 ) stop('SOME BIG ERROR'); summary(object) }}
res <- nmf(x, 3, nrun=3, .callback=cb())

## PARALLEL COMPUTATIONS
# try using 3 cores, but use sequential if not possible
res <- nmf(x, 3, nrun=3, .options='p3')

# force using 3 cores, error if not possible
res <- nmf(x, 3, nrun=3, .options='P3')

# use externally defined cluster
library(parallel)
cl <- makeCluster(6)
res <- nmf(x, 3, nrun=3, .pbackend=cl)

# use externally registered backend
registerDoParallel(cl)
res <- nmf(x, 3, nrun=3, .pbackend=NULL)
























#' # random data with underlying NMF model
v <- syntheticNMF(20, 3, 10)
# estimate a model
x <- nmf(v, 3)

# highligh row only (using custom colors)
basismap(x, tracks=':basis', annColor=list(basis=1:3))

## character annotation vector: ok if it does not contain 'basis'
# annotate first and second row + automatic special track
basismap(x, annRow=c('alpha', 'beta'))
# no special track here
basismap(x, annRow=c('alpha', 'beta', ':basis'), tracks=NA)
# with special track `basis`
basismap(x, annRow=list(c('alpha', 'beta'), ':basis'), tracks=NA)
# highligh columns only (using custom colors)
basismap(x, tracks='basis:')

# changing the name of the basis annotation track
basismap(x, annRow=list(new_name=':basis'))

# coefficient matrix
coefmap(x, annCol=c('alpha', 'beta')) # annotate first and second sample
coefmap(x, annCol=list('basis', Greek=c('alpha', 'beta'))) # annotate first and second sample + basis annotation
coefmap(x, annCol=c(new_name='basis'))





















# Generate random data
n <- 50; p <- 20
x <- abs(rmatrix(n, p, rnorm, mean=4, sd=1))
x[1:10, seq(1, 10, 2)] <- x[1:10, seq(1, 10, 2)] + 3
x[11:20, seq(2, 10, 2)] <- x[11:20, seq(2, 10, 2)] + 2
rownames(x) <- paste("ROW", 1:n)
colnames(x) <- paste("COL", 1:p)

## Scaling
aheatmap(x, scale = "row")
aheatmap(x, scale = "col") # partially matched to 'column'
aheatmap(x, scale = "r1") # each row sum up to 1
aheatmap(x, scale = "c1") # each colum sum up to 1

## Heatmap colors
aheatmap(x, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# color specification as an integer: use R basic colors
aheatmap(x, color = 1L)
# color specification as a negative integer: use reverse basic palette
aheatmap(x, color = -1L)
# color specification as a numeric: use HCL color
aheatmap(x, color = 1)
# do not cluster the rows
aheatmap(x, Rowv = NA)
# no heatmap legend
aheatmap(x, legend = FALSE)
# cell and font size
aheatmap(x, cellwidth = 10, cellheight = 5)

# directly write into a file
aheatmap(x, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "aheatmap.pdf")
unlink('aheatmap.pdf')

# Generate column annotations
annotation = data.frame(Var1 = factor(1:p %% 2 == 0, labels = c("Class1", "Class2")), Var2 = 1:10)

aheatmap(x, annCol = annotation)
aheatmap(x, annCol = annotation, annLegend = FALSE)


# Specify colors
Var1 = c("navy", "darkgreen")
names(Var1) = c("Class1", "Class2")
Var2 = c("lightgreen", "navy")

ann_colors = list(Var1 = Var1, Var2 = Var2)

aheatmap(x, annCol = annotation, annColors = ann_colors)

# Specifying clustering from distance matrix
drows = dist(x, method = "minkowski")
dcols = dist(t(x), method = "minkowski")
aheatmap(x, Rowv = drows, Colv = dcols)

# Display text in each cells
t <- outer(as.character(outer(letters, letters, paste0)), letters, paste0)[1:n, 1:p]
aheatmap(x, txt = t)
# NA values are shown as empty cells
t.na <- t
t.na[sample(length(t.na), 500)] <- NA # half of the cells
aheatmap(x, txt = t.na)

