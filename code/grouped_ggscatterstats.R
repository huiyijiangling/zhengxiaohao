rm(list=ls())
options(stringsAsFactors = F)
gc()
set.seed(123)
library(ggstatsplot)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)

load(r"(C:\Users\zxh\Desktop\R\meta collect\estimate_net25.Rdata)")
# basic function call
# ENSG00000140718 FTO
# ENSG00000091542 alkbh5
pdf("f2sfto_p.pdf",height=35,width = 20)
grouped_ggscatterstats(
  data             = estimate_net,
  x                = `ENSG00000140718`,
  y                = Purity,
  type             = "robust",
  plotgrid.args = list(ncol = 4),
  grouping.var     = project,
  ggplot.component = list(geom_rug(sides = "b")),
  annotation.args = list(tag_levels = "A")
)
dev.off()

pdf("f2sfto_sto.pdf",height=35,width = 20)

grouped_ggscatterstats(
  data             = estimate_net,
  x                = `ENSG00000140718`,
  y                = Stromal_score,
  type             = "robust",
  plotgrid.args = list(ncol = 4),
  grouping.var     = project,
  ggplot.component = list(geom_rug(sides = "b")),
  annotation.args = list(tag_levels = "A")
)
dev.off()

pdf("f2sfto_imm.pdf",height=35,width = 20)

grouped_ggscatterstats(
  data             = estimate_net,
  x                = `ENSG00000140718`,
  y                = Immune_score,
  type             = "robust",
  plotgrid.args = list(ncol = 4),
  grouping.var     = project,
  ggplot.component = list(geom_rug(sides = "b")),
  annotation.args = list(tag_levels = "A")
)
dev.off()

##########################


pdf("f2salkbh5_p.pdf",height=35,width = 20)
grouped_ggscatterstats(
  data             = estimate_net,
  x                = `ENSG00000091542`,
  y                = Purity,
  type             = "robust",
  plotgrid.args = list(ncol = 4),
  grouping.var     = project,
  ggplot.component = list(geom_rug(sides = "b")),
  annotation.args = list(tag_levels = "A")
)
dev.off()

pdf("f2salkbh5_sto.pdf",height=35,width = 20)
grouped_ggscatterstats(
  data             = estimate_net,
  x                = `ENSG00000091542`,
  y                = Stromal_score,
  type             = "robust",
  plotgrid.args = list(ncol = 4),
  grouping.var     = project,
  ggplot.component = list(geom_rug(sides = "b")),
  annotation.args = list(tag_levels = "A")
)
dev.off()

pdf("f2salkbh5_imm.pdf",height=35,width = 20)
grouped_ggscatterstats(
  data             = estimate_net,
  x                = `ENSG00000091542`,
  y                = Immune_score,
  type             = "robust",
  plotgrid.args = list(ncol = 4),
  grouping.var     = project,
  ggplot.component = list(geom_rug(sides = "b")),
  annotation.args = list(tag_levels = "A")
)
dev.off()

14*8