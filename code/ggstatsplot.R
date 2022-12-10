set.seed(123)
library(ggplot2)
library(ggstatsplot)
# https://zhuanlan.zhihu.com/p/107033510
# https://statsandr.com/blog/one-proportion-and-goodness-of-fit-test-in-r-and-by-hand/#assumption-of-prop.test-and-binom.test
## plot
#0
grouped_ggbarstats(
  data             = movies_long,
  x                = our,
  y                = mpaa,
  grouping.var  = cut,
  title            = "MPAA Ratings by Genre",
  xlab             = "movie genre",
  legend.title     = "MPAA rating",
  ggtheme          = hrbrthemes::theme_ipsum_pub(),
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  palette          = "Set2"
)
#1
ggbarstats(
  data             = movies_long,
  x                = genre,
  y                = mpaa,
  title            = "MPAA Ratings by Genre",
  xlab             = "movie genre",
  legend.title     = "MPAA rating",
  ggtheme          = hrbrthemes::theme_ipsum_pub(),
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  palette          = "Set2"
)
#2
ggbarstats(
  data         = mtcars,
  x            = am,
  y            = cyl,
  title            = "MPAA Ratings by Genre",
  xlab             = "movie genre",
  legend.title     = "MPAA rating",
  ggtheme          = hrbrthemes::theme_ipsum_pub(),
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  palette          = "Set2"
)
#3
ggpiestats(
  data         = mtcars,
  x            = am,
  y            = cyl,
  package      = "wesanderson",
  palette      = "Royal1",
  title        = "Dataset: Motor Trend Car Road Tests", ## title for the plot
  legend.title = "Transmission", ## title for the legend
  caption      = "Source: 1974 Motor Trend US magazine"
)
#4
ggpiestats(
  data             = movies_long,
  x                = genre,
  y                = mpaa,
  package      = "wesanderson",
  palette      = "Royal1",
  title        = "Dataset: Motor Trend Car Road Tests", ## title for the plot
  legend.title = "Transmission", ## title for the legend
  caption      = "Source: 1974 Motor Trend US magazine"
)
##
set.seed(123)

# plot
ggstatsplot::grouped_ggpiestats(
  data = ggstatsplot::movies_long,
  x = genre,
  grouping.var = mpaa, # grouping variable
  title.prefix = "Movie genre", # prefix for the faceted title
  label.repel = TRUE, # repel labels (helpful for overlapping labels)
  package = "ggthemr",
  palette = "dust",
  title.text = "Composition of MPAA ratings for different genres"
)

##
fisher.test(table(dat$Species, dat$size))
# 
# ggbarstats 这个更准确 但是直接上来的就是 
T1 <- xtabs(~genre + mpaa,data = movies_long)
T1
chisq.test(T1)$statistic
eeee$parameter[[1]]#DF
chisq.test(T1)$p.value


test <- fisher.test(table(dat$Species, dat$size))
test$p.value
test <- chisq.test(table(dat$Species, dat$size))
test$p.value

##################
test <- prop.test(
  x = 67, # number of heads
  n = 100, # number of trials
  p = 0.5 # expected probability of heads
)

test
##################
test <- binom.test(
  x = 12, # counts of successes
  n = 15, # total counts (12 + 3)
  p = 0.5 # expected proportion
)

test