
m6a_genelist=readxl::read_excel(r"(C:\Users\zxh\Desktop\杜永星老师 提交\m6a.xlsx)")
library(monocle3)
# cds_3d <- reduce_dimension(cds, max_components = 3)
# cds_3d <- cluster_cells(cds_3d)
# cds_3d <- learn_graph(cds_3d)
# cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
# cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")
# plot track
# plot_pseudotime_heatmap
scRNA=sce
data <- GetAssayData(scRNA, assay = 'RNA', slot = 'counts')
cell_metadata <- scRNA@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data),row.names=rownames(data))
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30,method="PCA")
# cds <- align_cds(cds)#,alignment_group = "batch")
cds <- reduce_dimension(cds,reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
# get_earliest_principal_node <- function(cds, time_bin="130-170"){
#   cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
#   
#   closest_vertex <-
#     cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <-
#     igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
#                                                               (which.max(table(closest_vertex[cell_ids,]))))]
#   
#   root_pr_nodes
# }
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
cds <- order_cells(cds)
# plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")
# #
plot_cells(cds,
           color_cells_by = "seurat_clusters",#"cell.type",
           genes="TFF1",#m6a_genelist$gene,#ciliated_genes,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = "seurat_clusters",#"cell.type",
           genes="TFF3",#m6a_genelist$gene,#ciliated_genes,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = "seurat_clusters",#"cell.type",
           genes="TFF1",#m6a_genelist$gene,#ciliated_genes,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)


plot_cells(cds,
           color_cells_by = "seurat_clusters",#"cell.type",
           genes="KRT19",#m6a_genelist$gene,#ciliated_genes,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)


plot_cells(cds,
           # color_cells_by = "embryo.time.bin",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

AFD_genes=m6a_genelist$gene
# AFD_genes <- c("gcy-8", "dac-1", "oig-8")
# AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
#                        colData(cds)$cell.type %in% c("AFD")]
AFD_genes=m6a_genelist$gene
# AFD_genes="FTO"
# AFD_genes="COL11A1"
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,]
# The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:

save(cds,file = "cds.Rdata")

pdf("m6a_genelist track.pdf",height = 8,width = 20)
p555=plot_genes_in_pseudotime(AFD_lineage_cds,ncol = 8,
                         # color_cells_by="embryo.time.bin",
                         min_expr=0.001) 
p555
dev.off()

# order.genes <- order.genes[!grepl("ENSMPUG", order.genes)]
# pdf("Fig5e.MP_M1_pseudotime_heatmap.pdf", 5/2.54*2,9/2.54*1.5)
# monocle::plot_pseudotime_heatmap(cds[rownames(cds)%in%AFD_genes,], num_clusters = 3, cores = 1, show_rownames = F, return_heatmap = T)
# dev.off()
pdf("m6a_genelist track111.pdf",height = 8,width = 20)
p555=plot_genes_in_pseudotime(cds,ncol = 8,
                              # color_cells_by="embryo.time.bin",
                              min_expr=0.001) 
p555
dev.off()


library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)

modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

print(htkm)
print(hthc)

# https://github.com/cole-trapnell-lab/monocle-release/issues/295
