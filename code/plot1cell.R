
library(plot1cell)
iri.integrated <- Install.example() 
save(iri.integrated,file="iri.integrated.Rdata")
load("iri.integrated.Rdata")

###Check and see the meta data info on your Seurat object
colnames(iri.integrated@meta.data)  

###Prepare data for ploting
circ_data <- prepare_circlize_data(iri.integrated, scale = 0.8 )
set.seed(1234)
cluster_colors<-rand_color(length(levels(iri.integrated)))
group_colors<-rand_color(length(names(table(iri.integrated$Group))))
rep_colors<-rand_color(length(names(table(iri.integrated$orig.ident))))

###plot and save figures
png(filename =  'circlize_plot.png', width = 6, height = 6,units = 'in', res = 300)
plot_circlize(circ_data,do.label = T, pt.size = 0.01, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
add_track(circ_data, group = "Group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()