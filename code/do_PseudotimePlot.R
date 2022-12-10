sample <- sce
# Transform into CDS.
cds <- SeuratWrappers::as.cell_data_set(sample)
pseudotime_genes="FTO"
pseudotime_genes="ALKBH5"
out <- SCpubr::do_PseudotimePlot(sample = sample,
                                 cds = cds,
                                 compute_monocle_partitions = TRUE,
                                 compute_monocle_clusters = FALSE,
                                 pseudotime_genes = pseudotime_genes)
# Retrieve trajectory groups.
p1 <- out$trajectory_groups
# Retrieve trajectory partitions.
p2 <- out$trajectory_partitions

p <- p1 | p2
p
