library(infercnv)
?CreateInfercnvObject
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
                                    annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="try2", 
                             cluster_by_groups=F, 
                             tumor_subcluster_partition_method="random_trees",
                             analysis_mode="subclusters",
                             denoise=TRUE,
                             HMM=TRUE,
                             num_threads = 4)



DATA=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv")


data(infercnv_data_example)#counts
data(infercnv_annots_example)


# V2
# normal_1
# normal
# normal_2
# normal
# normal_3
# normal
# normal_4
# normal
# normal_5
# normal
# normal_6
# normal
# normal_7
# normal
# normal_8
# normal
# normal_9
# normal
# normal_10
# normal
# tumor_grp_1_cell_1
# tumor
# tumor_grp_1_cell_2
# tumor
# tumor_grp_1_cell_3
# tumor
# tumor_grp_1_cell_4
# tumor
# tumor_grp_1_cell_5
# tumor
# tumor_grp_1_cell_6
# tumor
# tumor_grp_1_cell_7
# tumor
# tumor_grp_1_cell_8
# tumor
# tumor_grp_1_cell_9
# tumor
# tumor_grp_1_cell_10
# tumor

data(infercnv_genes_example)
# WASH7P chr1 14363 29806 anno