########## 生成 HIGH LOW 生成分组
# 10 Group-wise DE analysis plots
# 11 Term Enrichment Plots
# 15 Ligand-Receptor analysis
# 17 Pathway Activity inference analysis
# 18 TF Activity inference analysis
1 Dim plots
2 Feature plots

sce$assignment <- ifelse(sce$seurat_clusters %in% c("0", "2", "4"), "A", "B")

sce$assignment <- ifelse(sce$seurat_clusters %in% c("0", "2", "4"), "A", "B")

# Modify the space between nodes.
p <- SCpubr::do_SankeyPlot(sample = sce,
                           first_group = "seurat_clusters", #sub.cluster
                           middle_groups = c("assignment"),
                           last_group = "orig.ident",
                           type = "sankey",
                           colors.first = stats::setNames(
                             SCpubr::do_ColorPalette(colors.use = "steelblue",
                                                     n = length(unique(sce$sub.cluster))),
                             unique(sce$sub.cluster)
                           ))
p



p <- SCpubr::do_FeaturePlot(sample = sample,
                            features = c("nCount_RNA",
                                         "nFeature_RNA",
                                         "percent.mt",
                                         "CD14"),
                            plot.title = "A collection of features",
                            individual.titles = c("Plot A",
                                                  "Plot_B",
                                                  NA,
                                                  "Plot_D"),
                            ncol = 2)

p


3 Nebulosa plots
3.2 Compute joint densities

p <- SCpubr::do_NebulosaPlot(sample = sample, 
                             features = c("CD14", "CD8A"), 
                             joint = TRUE, 
                             return_only_joint = TRUE,
                             # individual.titles = c("Plot A",
                             #                       NA,
                             #                       "Combined density"),
                             plot.title = "Joint density CD14+-CD8A+")
p

10 Group-wise DE analysis plots
#################
# Seurat sample.
sample <- your_seurat_object

# Set the identities correctly.
Seurat::Idents(sample) <- sample$seurat_clusters

# Compute DE genes and transform to a tibble.
de_genes <- tibble::tibble(Seurat::FindAllMarkers(object = sample))
This is the basic output of the function:

# Default output.
p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                de_genes = de_genes)
##################
11 Term Enrichment Plots
genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")
p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                   dbs_use = "GO_Biological_Process_2021",
                                   nterms = 15)

12 Enrichment score heatmaps
A very common approach to make sense of your cells is to query several list of marker genes, retrieved from literature, and compute how enriched each cell is in each given list of genes. This is achieved by using Seurat::AddModuleScore. The scores can be then visualized as a Feature plot, but one can also aggregate the enrichment scores by any variable of interest, for instance the different clusters in the sample.
This kind of heatmaps can be easily computed using SCpubr::do_EnrichmentHeatmap():

genes <- list("Naive CD4+ T" = c("IL7R", "CCR7"),
              "CD14+ Mono" = c("CD14", "LYZ"),
              "Memory CD4+" = c("S100A4"),
              "B" = c("MS4A1"),
              "CD8+ T" = c("CD8A"),
              "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
              "NK" = c("GNLY", "NKG7"),
              "DC" = c("FCER1A", "CST3"),
              "Platelet" = c("PPBP"))

# Default parameters.
p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                  list_genes = genes)
p
###########
15 Ligand-Receptor analysis
One very interesting analysis that can be carried out is Ligand-Receptor analysis. This allows to compute whether specific clusters interact with each other based on the co-expression of a ligand and its receptor in the respective clusters. The interactions are retrieved from different databases and a plethora of tools have been released to tackle this analysis. One of them is liana, which is a framework that allows to run and integrate the results of several tools, providing a meta-analysis of the co-expression of ligand-receptor pairs. SCpubr makes use of liana and has its analysis and visualization integrated in the SCpubr::do_LigandReceptorPlot() function.

By default, the user has to run liana on their own and provide the resulting output as input for the function. The following code would produces the object that SCpubr::do_LigandReceptorPlot() expects as input:

liana_output <- liana::liana_wrap(sce = sample,
                                  method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
                                  idents_col = NULL,
                                  verbose = FALSE,
                                  assay = "SCT")
It is very important to note that liana_output has to contain the five different methods. This is a design choice. One can get the full output from scratch also running:

p <- SCPubr::do_LigandReceptorPlot(sample = sample)
The same output can be retrieved from:

# Ligand Receptor analysis plot.
p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output)
p
SCpubr do_LigandReceptorPlot default output.
Figure 1.3: SCpubr do_LigandReceptorPlot default output.

By default, top 25 unique, most significant interactions are retrieved and plotted. However, this can be changed by using top_interactions. Also, clusters that have no interactions, both as source and target, will be removed.:

# Ligand Receptor analysis plot with extra interactions.
p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                   top_interactions = 150)

######
16 Copy Number Variant analysis plots
Copy Number Variant analysis are another of the common analysis that one can compute on single-cell transcriptomics data. Provided with a reference one can compute, for the rest of the cells, whether there are Copy Number Variations (CNVs) in the cells across the chromosomes. This comes really handy to distinguish between tumor and healthy cells, provided that one has a CNV reference event to rely on. There are many tools to compute such analysis, but one that is highly used and that will be covered in this section is inferCNV.

InferCNV returns many outputs. The most interesting and straightforward one is the final image, such as the one below:

In it, we can observe that, for the different clusters, regions called as chromosome gains are colored as red and regions called as chromosome losses are colored as blue. This image, in the end, is a heatmap, meaning that there is, in the output of inferCNV, a matrix with the scores for each cell and gene, that we can make use of. This object is called, by default, run.final.infercnv_obj. This is the one that we will use.

16.1 Transferring the scores to a FeaturePlot.
One of the cool things we can do with this object, is to transfer the inferCNV scores back to our Seurat object and then plot them as a FeaturePlot. This can be achieved with the function SCpubr::do_CopyNumberVariantPlot(). For this function, one needs to provide the Seurat object and the final inferCNV object, together with the chromosome locations. If metacells were computed (not necessary, but used in this example), the mapping of cell-metacell has to be provided as well:

# This loads "human_chr_locations" into the environment.
utils::data("human_chr_locations", package = "SCpubr")

out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                        infercnv_object = infercnv_object,
                                        using_metacells = T,
                                        metacell_mapping = sample$metacell_mapping,
                                        chromosome_locations = human_chr_locations)
p <- out$`11p_umap`
p

##########
17 Pathway Activity inference analysis
This section is highly similar to TF Activity inference analysis as it makes use of the same source package and visualizations.

Another very important analysis that can be carried out, once we have our cells defined into groups (i.e: cell clusters), is Pathway Activity inference analysis. Out of the different tools out there to perform this task, the one we will use is called decoupleR (P.Badia-i-Mompel et al. 2022). This tool allows the inference of biological activities out of Omics data. For this, it requires a dataset to act as prior knowledge. For Pathway Activity inference, progeny is used. This allows for the computation of activity scores on a cell basis depicting how much (or little) each cell is enriched in each of the cancer Pathways stored in the database.

In order to visualize the enrichment of our cells in the pathways, the results need to be computed:

# Define your sample and assay.
sample <- your_seurat_object
assay <- "your_normalized_data_assay"

# Retrieve prior knowledge network.
network <- decoupleR::get_progeny(organism = "human")

# Run weighted means algorithm.
activities <- decoupleR::run_wmean(mat = as.matrix(sample@assays[[assay]]@data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "weight",
                                   times = 100,
                                   minsize = 5)


#####
18 TF Activity inference analysis
This section is highly similar to Pathway Activity inference analysis as it makes use of the same source package and visualizations.

Transcriptional Factor (TF) Activity inference analysis can also be carried out in a SC dataset. Same as with the previous chapter, we will use decoupleR (P.Badia-i-Mompel et al. 2022). For this, dorothea is used. This allows for the computation of activity scores on a cell basis depicting how much (or little) each cell is enriched in a TF and its dowstream targets (regulon).

In order to visualize the enrichment of our cells in the pathways, the results need to be computed:

# Define your sample and assay.
sample <- your_seurat_object
assay <- "your_normalized_data_assay"

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
With this, we can proceed to plot our results with SCpubr::do_TFActivityPlot().

18.1 Heatmap of averaged scores.

#######
19 Pseudotime analysis
Pseudotime analysis is, perhaps, one of the major analysis that can be carried out in SC data. It aims to depict a trajectory, a sort of order in which cells transition from A to B. What defines this process, and the starting and end point is heavily driven by the research question, the nature of the data and prior knowledge. While inferring this using data-driven approaches is possible, in SCpubr we will implement a prior-knowledge driven approach in a very specific use case.

For this, we will make use of monocle3. Out of the many uses of monocle3, we can infer a trajectory graph for the cells and, based on a list of marker genes, we can order the cells along this graph selecting either the cell with the highest or lowest enrichment score for the provided list of genes as starting point. For example:

In a general sense, we are working with a data set for which a given biological phenotype is expected. For this phenotype, we do have a list of genes that accurately characterize this process. The process needs to be of a differentiation kind, meaning that cells will transition based on this list of genes between state A and state B. Once these conditions are met, we can make use of the trajectory graph alongside enrichment scores for the cells in the list of genes. Depending whether this list of genes depict the start of the end of the process, we can use the cell with the highest or lowest enrichment scores as the starting point of the differentiation trajectory, ordering the rest of the cells based on that and computing the associated pseudotime.

Of course, selecting the one cell with the highest or lowest enrichment scores might be a simplistic approach, as this elevated or low enrichment score might be due to technical bias. In order to correct this a bit, we first select the cluster (or selected group) that, in average, has the highest or lowest mean enrichment score and from this cluster we select the cell with the highest or lowest score, accordingly.

Having this clear, again, we would like to remark that this is a very specific use case of this analysis and much more can be done using monocle3. Please have a look at their vignette for more information about pseudotime analysis. For our specific use case, SCpubr::do_PseudotimePlot() can be used.

19.1 Setting up partitions and clusters.
For this function to work, we need a Seurat object together with a Cell Data Set (CDS) object, that needs to be generated by SeuratWrappers::as.cell_data_set().

First thing to take into account is that monocle3 computes partitions and clusters. This is, when cells are too far away, the trajectory is disconnected, generating partitions with different start and end trajectories. We can tweak this process by forcing all of our data to be in a single partition or to keep as identities the group that we desire. This can be done by using compute_monocle_partitions and compute_monocle_clusters parameters.

We will use the following genes as a example, that depict CD14+ mono cells.

# Genes to use.
pseudotime_genes <- c("CD14", "LYN")
# Define your sample.
sample <- your_seurat_object
# Transform into CDS.
cds <- SeuratWrappers::as.cell_data_set(sample)

# Compute monocle clusters and partitions.
out <- SCpubr::do_PseudotimePlot(sample = sample,
                                 cds = cds,
                                 compute_monocle_partitions = TRUE,
                                 compute_monocle_clusters = TRUE,
                                 pseudotime_genes = pseudotime_genes)

# Compute monocle clusters and keep a single partition.
out <- SCpubr::do_PseudotimePlot(sample = sample,
                                 cds = cds,
                                 compute_monocle_partitions = FALSE,
                                 compute_monocle_clusters = TRUE,
                                 pseudotime_genes = pseudotime_genes)
#######
Save the figures
Creating good plots is just half of the process. It is equally important to properly save them. This is the purpose of SCpubr::save_Plot. This function is a very handy tool to save your plots easily in different formats, such as .pdf, .png, .jpeg, .tiff and .svg. This can be achieved by providing the following to output_format parameter:

"all": This will store the provided plot in all 5 formats.
"publication": This will store the plot in .pdf, .png and .svg.
Individual format: Provide the desired format and it will only be saved on that one.
Width and Height are set by default to 8 inches each, so the plot is squared. However, it is really important that these parameters are modified to the user’s need. The name of the file can be provided with file_name parameter and the path to store the files can be specified in figure_path. If not provided, figure_path will default to the current working environment and file_name will default to a combination of the current date and time. Here are some examples.

# Generate a plot.
p <- SCpubr::do_DimPlot(sample = sample)

# Default parameters.
SCpubr::save_Plot(plot = p)

# Specifying the name and folder.
SCpubr::save_Plot(plot = p,
                 figure_path = "/path/to/my/figures/",
                 file_name = "my_figure")

# Specify to also create a new folder.
SCpubr::save_Plot(plot = p,
                 figure_path = "/path/to/my/figures/",
                 file_name = "my_figure",
                 create_path = TRUE)

# Set dimensions for the figure.
SCpubr::save_Plot(plot = p,
                 figure_path = "/path/to/my/figures/",
                 file_name = "my_figure",
                 create_path = TRUE,
                 width = 8,
                 height = 8)

# Set quality (dpi).
SCpubr::save_Plot(plot = p,
                 figure_path = "/path/to/my/figures/",
                 file_name = "my_figure",
                 create_path = TRUE,
                 width = 8,
                 height = 8,
                 dpi = 300)





                 