
library(circlize)
set.seed(999)

bed = generateRandomBed()
head(bed)
bed = generateRandomBed(nr = 200, nc = 4)#nc 值
head(bed)
nrow(bed)
bed = generateRandomBed(nc = 2, fun = function(k) sample(letters, k, replace = TRUE))
head(bed)
####
circos.initializeWithIdeogram(species = "hg38")
text(0, 0, "default", cex = 1)
circos.info()
circos.clear()
circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(3,5,2,8)))
text(0, 0, "subset of chromosomes", cex = 1)



circos.par("start.degree" = 90)
circos.initializeWithIdeogram()
circos.clear()
text(0, 0, "'start.degree' = 90", cex = 1)
#######################
df = data.frame(
  name  = c("TP53",  "TP63",    "TP73"),
  start = c(7565097, 189349205, 3569084),
  end   = c(7590856, 189615068, 3652765))
circos.genomicInitialize(df)
tp_family = readRDS(system.file(package = "circlize", "extdata", "tp_family_df.rds"))
head(tp_family)
circos.genomicInitialize(tp_family)
circos.track(ylim = c(0, 1), 
             bg.col = c("#FF000040", "#00FF0040", "#0000FF40"), 
             bg.border = NA, track.height = 0.05)
n = max(tapply(tp_family$transcript, tp_family$gene, function(x) length(unique(x))))
# circos.genomicTrack(tp_family, ylim = c(0.5, n + 0.5), 
#                     panel.fun = function(region, value, ...) {
#                       all_tx = unique(value$transcript)
#                       for(i in seq_along(all_tx)) {
#                         l = value$transcript == all_tx[i]
#                         # for each transcript
#                         current_tx_start = min(region[l, 1])
#                         current_tx_end = max(region[l, 2])
#                         circos.lines(c(current_tx_start, current_tx_end), 
#                                      c(n - i + 1, n - i + 1), col = "#CCCCCC")
#                         circos.genomicRect(region[l, , drop = FALSE], ytop = n - i + 1 + 0.4, 
#                                            ybottom = n - i + 1 - 0.4, col = "orange", border = NA)
#                       }
#                     }, bg.border = NA, track.height = 0.4)
# circos.clear()
#############################
bed = generateRandomBed(nc = 2)
head(bed, n = 2)
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  if(CELL_META$sector.index == "chr1") {
    print(head(region, n = 2))
    print(head(value, n = 2))
  }
})
# circos.genomicTrackPlotRegion(data, ylim = c(0, 1),
#                               panel.fun = function(region, value, ...) {
#                                 circos.genomicPoints(region, value, ...)
#                               })
# circos.genomicTrackPlotRegion(data, numeric.column = c("value1", "value2"), 
#                               panel.fun = function(region, value, ...) {
#                                 circos.genomicPoints(region, value, ...)
#                               })
# circos.genomicPoints(region, value, numeric.column = c(1, 2))
# circos.genomicPoints(region, value, cex, pch)
# circos.genomicPoints(region, value, sector.index, track.index)
# circos.genomicTrack(data, numeric.column = 4, 
#                     panel.fun = function(region, value, ...) {
#                       # numeric.column is automatically passed to `circos.genomicPoints()`
#                       circos.genomicPoints(region, value, ...)
#                     })
# circos.genomicPoints = function(region, value, numeric.column = 1, ...) {
#   x = (region[[2]] + region[[1]])/2
#   for(i in numeric.column) {
#     y = value[[i]]
#     circos.points(x, y, ...)
#   }
# }
# circos.genomicLines(region, value, ...)
# circos.genomicLines(region, value, numeric.column = c(1, 2))
# circos.genomicLines(region, value, area, baseline, border)
# circos.genomicLines(region, value, sector.index, track.index)
# circos.genomicLines(region, value, lwd, lty = "segment")
# circos.genomicText(region, value, ...)
# circos.genomicText(region, value, y = 1, labels)
# circos.genomicText(region, value, numeric.column, labels.column)
# circos.genomicText(region, value, facing, niceFacing, adj)
# circos.genomicText(region, value, sector.index, track.index)
# circos.genomicRect(region, value, ytop = 1, ybottom = 0)
# circos.genomicRect(region, value, ytop.column = 2, ybottom = 0)
# circos.genomicRect(region, value, col, border)

set.seed(123)
bed1 = generateRandomBed(nr = 100)
bed1 = bed1[sample(nrow(bed1), 20), ]
bed2 = generateRandomBed(nr = 100)
bed2 = bed2[sample(nrow(bed2), 20), ]
?circos.genomicLink
circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), 
                   border = NA)
circos.clear()
# circos.genomicTrack(bed, ylim = c(-1, 1),
#                     panel.fun = function(region, value, ...) {
#                       circos.genomicPoints(region, value, ...)
#                       
#                       for(h in c(-1, -0.5, 0, 0.5, 1)) {
#                         circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
#                       }
#                       circos.text(x, y, labels)
#                       circos.axis("top")
#                     })
###############################
circos.initializeWithIdeogram(plotType = c("labels", "axis"))
circos.track(ylim = c(0, 1))
circos.genomicIdeogram() # put ideogram as the third track
circos.genomicIdeogram(track.height = 0.2)
#######################################################
circos.initializeWithIdeogram()
bed = generateRandomBed(nr = 100, nc = 4)
col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
circos.genomicHeatmap(bed, col = col_fun, side = "inside", border = "white")
circos.clear()
# circos.initializeWithIdeogram(plotType = NULL)
circos.initializeWithIdeogram()
circos.genomicHeatmap(bed, col = col_fun, side = "outside",
                      line_col = as.numeric(factor(bed[[1]])))
circos.genomicIdeogram()
circos.clear()
##############

circos.initializeWithIdeogram(,species = "hg38")
bed = generateRandomBed(nr = 50, nc = 4)
bed$symbolname=sample(letters, nrow(bed), replace = TRUE)
bed[1, 8] = "aaaaa"
circos.initializeWithIdeogram(plotType = c( "axis", "labels"))
circos.genomicLabels(bed, labels.column = ncol(bed), side = "inside",col = as.numeric(factor(bed[[1]])), line_col = as.numeric(factor(bed[[1]])))
col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
circos.genomicHeatmap(bed, numeric.column=4:7,col = col_fun, side = "inside", border = "white",connection_height=NULL,line_col = as.numeric(factor(bed[[1]])))
circos.link(bed1, 0, bed2, 0)
circos.clear()

set.seed(123)
bed1 = generateRandomBed(nr = 50,species = "hg38")# nc = 4,
set.seed(124)
# bed1 = bed1[sample(nrow(bed1), 20), ]
bed2 = generateRandomBed(nr = 50,species = "hg38")
set.seed(125)
# bed2 = bed2[sample(nrow(bed2), 20), ]
bed3 = generateRandomBed(nr = 50,species = "hg38")
bed3=dplyr::arrange(bed3,desc(chr))
# bed3 = bed3[sample(nrow(bed3), 20), ]
bed_list=list(bed1,bed2,bed3)
bed1$chr="chr1"
bed1$start=as.integer("12716970")
bed1$end=as.integer("83063822")
#F1 不标记名字展示全景
pdf("circos.pdf",height=8,width=8)
circos.initializeWithIdeogram(,species = "hg38")
circos.genomicRainfall(bed_list, col = c("red", "yellow","blue"),pch = 16, cex = 0.4,track.height = 0.1)
circos.genomicDensity(bed_list[[1]], col = c("red"), track.height = 0.1)
circos.genomicDensity(bed_list[[2]], col = c("yellow"), track.height = 0.1)
circos.genomicDensity(bed_list[[3]], col = c("blue"), track.height = 0.1)
if(F){
aaaa=lapply(1:10^8, function(k) rand_color(1,luminosity = "light", transparency = 0.5))
aaaa_color=unique(unlist(aaaa))
save(aaaa_color,file="color from rand color.Rdata")
load("color from rand color.Rdata")
aaaa_color[1:nrow(bed)]
circos.initializeWithIdeogram(plotType = NULL)
}
circos.genomicLink(bed1, bed2, col = aaaa_color[1:nrow(bed1)], lty=3,lwd = 0.001,
                   border = NA,directional = -1,arr.col="red",arr.lty=0,arr.length=0.2)#col = rand_color(nrow(bed1), transparency = 0.5)
circos.genomicLink(bed2, bed3, col =aaaa_color[1:nrow(bed1)],lty=3,lwd = 0.001, 
                   border = NA,directional = 1,arr.col="blue",arr.lty=0,arr.length=0.2)
dev.off()
circos.initializeWithIdeogram()
circos.genomicLabels(bed, labels.column = 4, side = "outside",
                     col = as.numeric(factor(bed[[1]])), line_col = as.numeric(factor(bed[[1]])))
circos.genomicIdeogram()
circos.clear()


circos.initializeWithIdeogram(plotType = c( "axis", "labels"))
circos.genomicIdeogram()
circos.axis()
circos.labels()
?circos.axis
# still work on the ideogram track
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "top")
})
circos.track(ylim = c(0, 1), track.height = 0.1)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "bottom", direction = "inside")
})
#
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicIdeogram()
# still work on the ideogram track
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "top")
})
circos.track(ylim = c(0, 1), track.height = 0.1)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "bottom", direction = "inside")
})
# circos.clear()
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicIdeogram()
# still work on the ideogram track
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "top")
})
circos.genomicRainfall(bed)
circos.genomicDensity(bed, baseline = 0)
circos.genomicDensity(bed, window.size = 1e6)
circos.genomicDensity(bedlist, col = c("#FF000080", "#0000FF80"))

load(system.file(package = "circlize", "extdata", "DMR.RData"))
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))

bed_list = list(DMR_hyper, DMR_hypo)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))
circos.genomicDensity(DMR_hyper, col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(DMR_hypo, col = c("#0000FF80"), track.height = 0.1)
circos.clear()
head(rainfallTransform(DMR_hyper))
head(genomicDensity(DMR_hyper, window.size = 1e6))
###########






load(system.file(package = "circlize", "extdata", "DMR.RData"))
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))

bed_list = list(DMR_hyper, DMR_hypo)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))
circos.genomicDensity(DMR_hyper, col = c("#FF000080"), track.height = 0.1)#count_by = c("percent", "number"),
circos.genomicDensity(DMR_hypo, col = c("#0000FF80"), track.height = 0.1)#count_by = c("percent", "number"),


circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(DMR_hyper, col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(DMR_hyper, col = c("#FF000080"), count_by = "number", track.height = 0.1)

