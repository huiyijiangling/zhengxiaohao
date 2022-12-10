sample <- sce
assay <- "RNA"

# Retrieve prior knowledge network.
network <- decoupleR::get_dorothea(organism = "human",
                                   levels = c("A", "B", "C"))

# Run weighted means algorithm.
activities <- decoupleR::run_wmean(mat = as.matrix(sample@assays[[assay]]@data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "mor",
                                   times = 100,
                                   minsize = 5)
save(activities,file="activities.Rdata")
sample$split.me <- ifelse(sample$seurat_clusters %in% c("0", "3", "7"), "Group A","Group B")

out <- SCpubr::do_TFActivityPlot(sample = sample,
                                 activities = activities)
                                 # ,
                                 # split.by = "split.me")
p <- out$heatmaps$average_scores
p

