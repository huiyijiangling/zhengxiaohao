expr$
single.savePaths; sampleNames; savePath; combName;
authorName = NULL;
comb.method = "NormalMNN";
harmony.theta = NULL;
harmony.lambda = NULL;
harmony.sigma = 0.1;
vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent","S.Score","G2M.Score");#"CellCycle.score");
pc.use = 30;
resolution = 0.8;
clusterStashName = "comb.cluster";
show.features = NULL; bool.add.features = T;
bool.runDiffExpr = T;
n.markers = 5;
sample.colors = NULL;
species = "human";
genome = "hg19";
hg.mm.mix = F;
bool.runCellClassify = T;
ct.templates = NULL;
coor.names = c("tSNE_1", "tSNE_2");
bool.runMalignancy = T;
cnv.ref.data = NULL;
cnv.referAdjMat = NULL;
cutoff = 0.1;
p.value.cutoff = 0.5;
bool.intraTumor = T;
bool.runCellCycle = T;
bool.runStemness = T;
bool.runGeneSets = T;
geneSets = NULL;
geneSet.method = "average";
bool.runExprProgram = T;
nmf.rank = 50;
genReport = T


if(T){
  single.savePaths <- list.files(r"(J:\cra\results\)",full.names = T,pattern = "T|N")
  sampleNames <- list.files(r"(J:\cra\results\)",pattern = "T|N")    # The labels for all samples
  savePath <- "./results/comb"       # A path to save the results
  combName <- "comb"                 # A label of the combined samples
  authorName <- "Xiaobei"                # The author name to mark the report
  comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
}

single.savePaths = single.savePaths;
sampleNames = sampleNames;
savePath = savePath;
combName = combName;
authorName = authorName;
pc.use = 30;
resolution = 0.3;#人cra0.3
bool.runDiffExpr = F;
bool.runMalignancy = F;
species ="human";#"mouse";#
genome = "hg38";#"mm10",#"hg19",
coor.names = c("tSNE_1", "tSNE_2");
# bool.runCellCycle = F;
# bool.runGeneSets = F;
# coor.names = c("UMAP_1", "UMAP_2");
comb.method = comb.method

table(expr.list[[1]]$mito.percent>0.1)

table(c(cc.genes$s.genes,cc.genes$g2m.genes)%in%VariableFeatures(expr))
length(VariableFeatures(expr))
