library(circlize)
library(ggplot2)
bed = generateRandomBed(nr = 200)
head(bed)
#    chr    start      end     value1
# 1 chr1  6255701 16825950 -0.7966037
# 2 chr1 23545777 26681342 -1.4405458
# 3 chr1 41443165 46719497  0.1430059
# 4 chr1 52909447 54931298 -0.1015327
# 5 chr1 65273586 74868674 -0.5058086
# 6 chr1 75096106 79815973  0.16868
然后第一步就是先绘制基因组，也是一个函数就可以实现，这里的基因组选择的是 hg38 ：
pdf("ggg.pdf",width = 8,height = 8)
circos.initializeWithIdeogram(species = "hg38")


circos.genomicTrack(bed, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = 1, ...)
                    })
                      
                      circos.genomicTrack(bed, 
                                          panel.fun = function(region, value, ...) {
                                            circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                                               col = ifelse(value[[1]] > 0, "red", "green"), ...)
                                            circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                                          })
                      
                      circos.clear()
dev.off()

                      maf = read.table('./7.anotation/vep/case1_biorep_A_techrep_vep.maf', comment.char = "#",header = T,sep = '\t',quote = "")
                      choose.col = c("Chromosome","Start_Position","End_Position","Hugo_Symbol","t_depth","t_ref_count","t_alt_count")
                      somatic = maf[ ,choose.col]
                      somatic$vaf = somatic$t_alt_count/somatic$t_depth
                      somatic.bed = somatic[,c(1,2,3,8)]
                      CNV
                      seg = read.table('./8.cnv/gatk/segment/case1_biorep_A_techrep.cr.igv.seg', header=T)
                      seg$Segment_Mean[ seg$Segment_Mean < -1 ] = -1
                      
                      seg.bed = seg[,c(2,3,4,6)]
                      可视化
                      circos.initializeWithIdeogram(species = "hg38")
                      
                      circos.genomicTrack(somatic.bed, 
                                          panel.fun = function(region, value, ...) {
                                            circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = 1, ...)
                                          })
                      circos.genomicTrack(seg.bed, 
                                          panel.fun = function(region, value, ...) {
                                            circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                                               col = ifelse(value[[1]] > 0, "red", "blue"), ...)
                                            circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                                          })
                      
                      circos.clear()