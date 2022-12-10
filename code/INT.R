library(RNOmni)
# Draw from chi-1 distribution
y <- rchisq(n = 1e3, df = 1)
# Rank normalize
z <- RankNorm(y)
# Plot density of transformed measurement
plot(density(z))
box-cox 是另一个方法  最后需要去 back transform 来增加解释性
