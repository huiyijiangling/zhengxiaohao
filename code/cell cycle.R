



# filtered_seurat=sce
seurat_phase$CellCycle.score
seurat_phase$G2M.Score
seurat_phase$S.Score

seurat_phase <- NormalizeData(filtered_seurat)
seurat_phase <- FindVariableFeatures(seurat_phase,
                                     nfeatures = 2000, selection.method = "vst")
seurat_phase <- ScaleData(seurat_phase, vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent","CC.Difference"))#
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes)
seurat_phase$CC.Difference <- seurat_phase$S.Score - seurat_phase$G2M.Score

# Score cells for cell cycle
# normalize之后就可以算了实际上教程是scale后
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes)
RidgePlot(seurat_phase, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
#normalize 无论后续怎么scale or vars.to.regress
# 理论上是scale后来比
# seurat_phase$S.Score["N1_N1_AAACCTGCAACCGCCA"]
# # 0.0939852 
# seurat_phase$G2M.Score["N1_N1_AAACCTGCAACCGCCA"]
# # 0.01974254
# seurat_phase$Phase["N1_N1_AAACCTGCAACCGCCA"]
# "S"

#ori
# seurat_phase1$S.Score["N1_N1_AAACCTGCAACCGCCA"]
# # 0.2023756 
# seurat_phase1$G2M.Score["N1_N1_AAACCTGCAACCGCCA"]
# # -0.0490566 
# seurat_phase1$Phase["N1_N1_AAACCTGCAACCGCCA"]
# "S"
# marrow$CC.Difference <- marrow$S.Score - marrow$G2M.Score

marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
seurat_phase <- RunPCA(seurat_phase, features = c(cc.genes$g2m.genes, cc.genes$s.genes))
DimPlot(seurat_phase,reduction = "pca")
