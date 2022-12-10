mmend using sctransform (Hafemeister and Satija, Genome Biology 2019), which which builds regularized negative binomial models of gene expression in order to account for technical artifacts while preserving biological variance. For more details on sctransform, please see the paper here and the Seurat vignette here. sctransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
How do results compare to log-normalization?
  Gene expression visualization
In Seurat, we have functionality to explore and interact with the inherently visual nature of spatial data. The SpatialFeaturePlot() function in Seurat extends FeaturePlot(), and can overlay molecular data on top of tissue histology. For example, in this data set of the mouse brain, the gene Hpca is a strong hippocampus marker and Ttr is a marker of the choroid plexus.

SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))


The default parameters in Seurat emphasize the visualization of molecular data. However, you can also adjust the size of the spots (and their transparency) to improve the visualization of the histology image, by changing the following parameters:
  
  pt.size.factor- This will scale the size of the spots. Default is 1.6
alpha - minimum and maximum transparency. Default is c(1, 1).
Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2


Dimensionality reduction, clustering, and visualization
We can then proceed to run dimensionality reduction and clustering on the RNA expression data, using the same workflow as we use for scRNA-seq analysis.

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
We can then visualize the results of the clustering either in UMAP space (with DimPlot()) or overlaid on the image with SpatialDimPlot().

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2
As there are many colors, it can be challenging to visualize which voxel belongs to which cluster. We have a few strategies to help with this. Setting the label parameter places a colored box at the median of each cluster (see the plot above).

You can also use the cells.highlight parameter to demarcate particular cells of interest on a SpatialDimPlot(). This can be very useful for distinguishing the spatial localization of individual clusters, as we show below:
  
  SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,
                                                                                       5, 8)), facet.highlight = TRUE, ncol = 3)


Interactive plotting
We have also built in a number of interactive plotting capabilities. Both SpatialDimPlot() and SpatialFeaturePlot() now have an interactive parameter, that when set to TRUE, will open up the Rstudio viewer pane with an interactive Shiny plot. The example below demonstrates an interactive SpatialDimPlot() in which you can hover over spots and view the cell name and current identity class (analogous to the previous do.hover behavior).

SpatialDimPlot(brain, interactive = TRUE)

For SpatialFeaturePlot(), setting interactive to TRUE brings up an interactive pane in which you can adjust the transparency of the spots, the point size, as well as the Assay and feature being plotted. After exploring the data, selecting the done button will return the last active plot as a ggplot object.

SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)

The LinkedDimPlot() function links the UMAP representation to the tissue image representation and allows for interactive selection. For example, you can select a region in the UMAP plot and the corresponding spots in the image representation will be highlighted.

LinkedDimPlot(brain)

Identification of Spatially Variable Features
Seurat offers two workflows to identify molecular features that correlate with spatial location within a tissue. The first is to perform differential expression based on pre-annotated anatomical regions within the tissue, which may be determined either from unsupervised clustering or prior knowledge. This strategy works will in this case, as the clusters above exhibit clear spatial restriction.

de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)


An alternative approach, implemented in FindSpatiallyVariables(), is to search for features exhibiting spatial patterning in the absence of pre-annotation. The default method (method = 'markvariogram), is inspired by the Trendsceek, which models spatial transcriptomics data as a mark point process and computes a ‘variogram’, which identifies genes whose expression level is dependent on their spatial location. More specifically, this process calculates gamma(r) values measuring the dependence between two spots a certain “r” distance apart. By default, we use an r-value of ‘5’ in these analyses, and only compute these values for variable genes (where variation is calculated independently of spatial location) to save time.

We note that there are multiple methods in the literature to accomplish this task, including SpatialDE, and Splotch. We encourage interested users to explore these methods, and hope to add support for them in the near future.

brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
    selection.method = "markvariogram")
Now we visualize the expression of the top 6 features identified by this measure.

top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))


Subset out anatomical regions
As with single-cell objects, you can subset the object to focus on a subset of data. Here, we approximately subset the frontal cortex. This process also facilitates the integration of these data with a cortical scRNA-seq dataset in the next section. First, we take a subset of clusters, and then further segment based on exact positions. After subsetting, we can visualize the cortical cells either on the full image, or a cropped image.

cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2


Integration with single-cell data
At ~50um, spots from the visium assay will encompass the expression profiles of multiple cells. For the growing list of systems where scRNA-seq data is available, users may be interested to ‘deconvolute’ each of the spatial voxels to predict the underlying composition of cell types. In preparing this vignette, we tested a wide variety of decovonlution and integration methods, using a reference scRNA-seq dataset of ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, generated with the SMART-Seq2 protocol. We consistently found superior performance using integration methods (as opposed to deconvolution methods), likely because of substantially different noise models that characterize spatial and single-cell datasets, and integration methods are specifiically designed to be robust to these differences. We therefore apply the ‘anchor’-based integration workflow introduced in Seurat v3, that enables the probabilistic transfer of annotations from a reference to a query set. We therefore follow the label transfer workflow introduced here, taking advantage of sctransform normalization, but anticipate new methods to be developed to accomplish this task.

We first load the data (download available here), pre-process the scRNA-seq reference, and then perform label transfer. The procedure outputs, for each spot, a probabilistic classification for each of the scRNA-seq derived classes. We add these predictions as a new assay in the Seurat object.

allen_reference <- readRDS("../data/allen_cortex.rds")
# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)


anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
    weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
Now we get prediction scores for each spot for each class. Of particular interest in the frontal cortex region are the laminar excitatory neurons. Here we can distinguish between distinct sequential layers of these neuronal subtypes, for example:

DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)


Based on these prediction scores, we can also predict cell types whose location is spatially restricted. We use the same methods based on marked point processes to define spatially variable features, but use the cell type prediction scores as the “marks” rather than gene expression.

cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "markvariogram",
    features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)


Finally, we show that our integrative procedure is capable of recovering the known spatial localization patterns of both neuronal and non-neuronal subsets, including laminar excitatory, layer-1 astrocytes, and the cortical grey matter.

SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT",
    "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))


Working with multiple slices in Seurat
This dataset of the mouse brain contains another slice corresponding to the other half of the brain. Here we read it in and perform the same initial normalization.

brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
In order to work with multiple slices in the same Seurat object, we provide the merge function.

brain.merge <- merge(brain, brain2)
This then enables joint dimensional reduction and clustering on the underlying RNA expression data.

DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)
Finally, the data can be jointly visualized in a single UMAP plot. SpatialDimPlot() and SpatialFeaturePlot() will by default plot all slices as columns and groupings/features as rows.

DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))


SpatialDimPlot(brain.merge)


SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))


Acknowledgments
We would like to thank Nigel Delaney and Stephen Williams for their helpful feedback and contributions to the new spatial Seurat code.

Slide-seq
Dataset
Here, we will be analyzing a dataset generated using Slide-seq v2 of the mouse hippocampus. This tutorial will follow much of the same structure as the spatial vignette for 10x Visium data but is tailored to give a demonstration specific to Slide-seq data.

You can use our SeuratData package for easy data access, as demonstrated below. After installing the dataset, you can type ?ssHippo to see the commands used to create the Seurat object.

InstallData("ssHippo")
slide.seq <- LoadData("ssHippo")
Data preprocessing
The initial preprocessing steps for the bead by gene expression data are similar to other spatial Seurat analyses and to typical scRNA-seq experiments. Here, we note that many beads contain particularly low UMI counts but choose to keep all detected beads for downstream analysis.

plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)


We then normalize the data using sctransform and perform a standard scRNA-seq dimensionality reduction and clustering workflow.

slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)
We can then visualize the results of the clustering either in UMAP space (with DimPlot()) or in the bead coordinate space with SpatialDimPlot().

plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2


SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, idents = c(1,
    6, 13)), facet.highlight = TRUE)


Integration with a scRNA-seq reference
To facilitate cell-type annotation of the Slide-seq dataset, we are leveraging an existing mouse single-cell RNA-seq hippocampus dataset, produced in Saunders*, Macosko*, et al. 2018. The data is available for download as a processed Seurat object here, with the raw count matrices available on the DropViz website.

ref <- readRDS("../data/mouse_hippocampus_reference.rds")
The original annotations from the paper are provided in the cell metadata of the Seurat object. These annotations are provided at several “resolutions”, from broad classes (ref$class) to subclusters within celltypes (ref$subcluster). For the purposes of this vignette, we’ll work off of a modification of the celltype annotations (ref$celltype) which we felt struck a good balance.

We’ll start by running the Seurat label transfer method to predict the major celltype for each bead.

anchors <- FindTransferAnchors(reference = ref, query = slide.seq, normalization.method = "SCT",
    npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype, prediction.assay = TRUE,
    weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay
We can then visualize the prediction scores for some of the major expected classes.

DefaultAssay(slide.seq) <- "predictions"
SpatialFeaturePlot(slide.seq, features = c("Dentate Principal cells", "CA3 Principal cells", "Entorhinal cortex",
    "Endothelial tip", "Ependymal", "Oligodendrocyte"), alpha = c(0.1, 1))


slide.seq$predicted.id <- GetTransferPredictions(slide.seq)
Idents(slide.seq) <- "predicted.id"
SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, idents = c("CA3 Principal cells",
    "Dentate Principal cells", "Endothelial tip")), facet.highlight = TRUE)


Identification of Spatially Variable Features
As mentioned in the Visium vignette, we can identify spatially variable features in two general ways: differential expression testing between pre-annotated anatomical regions or statistics that measure the dependence of a feature on spatial location.

Here, we demonstrate the latter with an implementation of Moran’s I available via FindSpatiallyVariableFeatures() by setting method = 'moransi'. Moran’s I computes an overall spatial autocorrelation and gives a statistic (similar to a correlation coefficient) that measures the dependence of a feature on spatial location. This allows us to rank features based on how spatially variable their expression is. In order to facilitate quick estimation of this statistic, we implemented a basic binning strategy that will draw a rectangular grid over Slide-seq puck and average the feature and location within each bin. The number of bins in the x and y direction are controlled by the x.cuts and y.cuts parameters respectively. Additionally, while not required, installing the optional Rfast2 package(install.packages('Rfast2')), will significantly decrease the runtime via a more efficient implementation.

DefaultAssay(slide.seq) <- "SCT"
slide.seq <- FindSpatiallyVariableFeatures(slide.seq, assay = "SCT", slot = "scale.data", features = VariableFeatures(slide.seq)[1:1000],
    selection.method = "moransi", x.cuts = 100, y.cuts = 100)
Now we visualize the expression of the top 6 features identified by Moran’s I.

SpatialFeaturePlot(slide.seq, features = head(SpatiallyVariableFeatures(slide.seq, selection.method = "moransi"),
    6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")