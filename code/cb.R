authorName = NULL;
comb.method = "NormalMNN";
harmony.theta = NULL;
harmony.lambda = NULL;
harmony.sigma = 0.1;
vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent");
pc.use = 30;
resolution = 0.8;
clusterStashName = "comb.cluster";
show.features = NULL; bool.add.features = T;
bool.runDiffExpr = F;
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


single.savePaths <- list.files(r"(H:\download.big.ac.cn\gsa\CRA001160\results\)",full.names = T,pattern = "T|N")
sampleNames <- list.files(r"(H:\download.big.ac.cn\gsa\CRA001160\results\)",pattern = "T|N")    # The labels for all samples
savePath <- "./results/comb"       # A path to save the results
combName <- "comb"                 # A label of the combined samples
authorName <- "Xiaobei"                # The author name to mark the report
comb.method <- "Harmony"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")