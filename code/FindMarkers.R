library(ComplexHeatmap)
sample <- sce
assay <- "RNA"

# Set the identities correctly.
Seurat::Idents(sce) <- sample$seurat_clusters

#######################

de_genes2=FindMarkers(sce_T,group.by = 'FTO_group',ident.1 = )
  # de_genes%>%
  # group_by(cluster) %>%
  # # mutate(var = var(avg)) %>%
  # # ungroup() %>%
  # top_n(5,avg_log2FC) %>%
  # distinct()
# Add more layers of mean expression with group.by.
p <- SCpubr::do_GroupwiseDEPlot(sample = sce,
                                de_genes = de_genes2)
                                ,
                                group.by = c("seurat_clusters","FTO_group"),
                                # group.by = c("seurat_clusters",
                                             # "modified_orig.ident", 
                                             # "orig.ident"),
                                row_title_expression = c("", # "Title A",
                                                         "Title B")
)

p
#########################
genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")
# “C”: Performs a query for the functional terms (MsigDB, GO-BP, GO-MF and KEGG).
p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                   dbs_use = "C")
patchwork::wrap_plots(p, nrow = 2)
