library(impute)
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load the package
library(WGCNA);
# Allow multi-threading in WGCNA
allowWGCNAThreads();
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load custom functions that simplify handling of multi-data sets
source("./RCode/CommonFunctions/networkFunctions-extras-05GAI.R");
# Base of the data tree
dataDir = "C:/Users/zxh/Desktop/R/wgcna_meta/Project-MetaAnalysis/LungCancer/Data-Expression/";
# Next level is split by data set
setDir = c("LungCancer-Shedden/MICH/",
           "LungCancer-Shedden/HLM/",
           "LungCancer-Shedden/DFCI/",
           "LungCancer-Shedden/MSKCC/",
           # "LungCancer-Bild-GSE03141/",
           "LungCancer-Tomida-GSE13213/",
           "LungCancer-Taekuchi-GSE11969/",
           "LungCancer-Roepman/")
# Under each data set, the relevant data set is stored in this subdirectory
exprDirNormalized = "Expression/020-normalized/"
# Actual data files
exprFiles = spaste(dataDir,
                   setDir,
                   exprDirNormalized,
                   c("GeneExpressionMICH.csv.bz2",
                     "GeneExpressionHLM.csv.bz2",
                     "GeneExpressionDFCI.csv.bz2",
                     "GeneExpressionMSKCC.csv.bz2",
                     "GeneExpressionGSE3141.csv.bz2",
                     "GeneExpressionTomida.csv.bz2",
                     "GeneExpressionTakeuchi.csv.bz2",
                     "GeneExpressionRoepman.csv.bz2"));
# Directory containing microarray platform annotation files
annotDir = spaste(dataDir, c("LungCancer-Shedden/", spaste(setDir[5:8], "Expression/")),
                  "ArrayAnnotation/Shortened/");
# Actual annotation files
annotFiles =spaste(dataDir, c("LungCancer-Shedden/", spaste(setDir[5:8], "Expression/")),
                   "ArrayAnnotation/",
                   c("Shortened/ShortGeneSymbolHG133A.csv.bz2",
                     "Shortened/ShortGeneSymbolHG133plus2.csv.bz2",
                     "Raw-Large/LongTomidaGPL6480-26599.txt.bz2",
                     "Raw-Large/GeneSymbolTakeuchi.csv.bz2",
                     "Shortened/Roepman_ArrayAnnotation.txt.bz2"))
# Index indicating which annotation files is to be used for each data set
annotInd = c(1,1,1,1,2,3,4,5);
nAnnot = length(annotFiles)
# Set names for pretty-printing
setNames = c("Shedden-MICH",
             "Shedden-HLM",
             "Shedden-DFCI",
             "Shedden-MSKCC",
             "Bild",
             "Tomida",
             "Takeuchi",
             "Roepman");

manAnnotFiles = spaste(annotDir,
                       c("HG-U133A.na31.annot-shortened.txt.bz2",
                         # "HG-U133_Plus_2.na31.annot-shortened.txt.bz2",
                         "Agilent014850_D_AA_20070207-Tomida-shortened.txt.bz2",
                         "AgilentCustomHGArray-Takeuchi.csv.bz2",
                         "Roepman_ArrayAnnotation.txt.bz2"));
type = c(1,1,1,2,1);
mannot = list();
for (a in 1:nAnnot)
{
  if (type[a]==2)
  {
    mannot[[a]] = read.csv(bzfile(manAnnotFiles[a]), header = TRUE, colClasses = "character");
  } else
    mannot[[a]] = read.delim(bzfile(manAnnotFiles[a]), header = TRUE, comment.char = "#",
                             colClasses = "character");
} 
#Quick look at what is included:
lapply(mannot, colnames)
lapply(mannot, dim)
idConversionFiles = spaste(annotDir,
                           c("Entrez-GBlist-HG133A-convertedByDavid.txt.bz2",
                             # "Entrez-GBlist-HG133plus2-convertedByDavid.txt.bz2",
                             "Entrez-GBlist-Tomida-convertedByDavid.txt.bz2",
                             "Entrez-GBlist-Takeuchi-ConvertedByDavid.txt.bz2",
                             "Entrez-HSlist-Roepman-convertedByDavid.txt.bz2"));
idConv = list();
# Read in the conversion files
for (a in 1:nAnnot)
{
  idConv[[a]] = read.delim(bzfile(idConversionFiles[a]), header = TRUE, colClasses = "character");
} 
#Remove the superfluous "0" identifiers from Roepman
a = 5
idConv[[a]] = idConv[[a]] [ idConv[[a]][, 1]!="0", ];
# We actually need only the last two id conversions.
mannot.ext = mannot;
a = 4;
table(is.finite(match(idConv[[a]]$From, mannot[[a]]$GB_LIST)))
mannot.ext[[a]] = merge(mannot[[a]], idConv[[a]], by.x = "GB_LIST", by.y = "From");
colnames(mannot.ext[[a]])[colnames(mannot.ext[[a]])=="To"] = "Entrez";
a = 5;
mannot[[a]]$GeneID = toupper(mannot[[a]]$GeneID);
idConv[[a]]$From = toupper(idConv[[a]]$From);
table(is.finite(match(idConv[[a]]$From, mannot[[a]]$GeneID)))