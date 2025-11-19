#
# Script: 02_run_WGCNA_Reproduce.R
#
# Description:
# Performs the full Weighted Gene Co-expression Network Analysis (WGCNA)
# following the exact logic of the user-provided original script
# for reproducibility.
#
# Project:      MM SRGs Prognostic Model
#
# ---

# ---
# Section 1: Load Libraries and Configuration
# ---
library(limma)
library(WGCNA)
options(stringsAsFactors = FALSE)

# Input Files
EXPRESSION_FILE <- "../data/input/wgcna_expression_input.txt" # Original: 工3SP代mrnasi与免疫打分.txt
CLINICAL_FILE   <- "../data/input/wgcna_clinical_traits.txt"  # Original: 工3打分列表.txt

# Output Directories
FIG_DIR     <- "../results/figures/wgcna/"
TABLE_DIR   <- "../results/tables/wgcna/"
MODULE_DIR  <- "../results/tables/wgcna/module_genes/"

# Create directories if they don't exist
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MODULE_DIR, showWarnings = FALSE, recursive = TRUE)

# Parameters (from original script)
CUT_HEIGHT <- 60
SOFT_POWER_BETA <- 11
MIN_MODULE_SIZE <- 100
MERGE_CUT_HEIGHT <- 0.05

# ---
# Section 2: Load Data
# ---
data <- read.table(EXPRESSION_FILE, header = TRUE, check.names = FALSE, row.names = 1)
data=data[rowMeans(data)>0.5,]

cli2=read.table(CLINICAL_FILE, header=T, sep="\t", check.names=F, row.names=1)

# ---
# Section 3: Preprocessing and Outlier Detection
# ---
datExpr0=t(data)

str(datExpr0) 

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average")

pdf(file = file.path(FIG_DIR, "01_sample_cluster_1.pdf"), width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

pdf(file = file.path(FIG_DIR, "01_sample_cluster_2.pdf"), width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = CUT_HEIGHT, col = "red")
dev.off()

# ---
# Section 4: Sample Filtering (Restored from original logic to prevent hclust error)(have been performed)
# ---
#clust = cutreeStatic(sampleTree, cutHeight = CUT_HEIGHT, minSize = 10)
#keepSamples = (clust==1)
#datExpr0 = datExpr0[keepSamples, ]

# ---
# Section 5: Soft Threshold Selection
# ---
enableWGCNAThreads()
powers = c(1:20)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

pdf(file=file.path(FIG_DIR, "02_scale_independence.pdf"),width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=cex1,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# ---
# Section 6: Network Construction and Module Detection
# ---
sft
power = SOFT_POWER_BETA
adjacency = adjacency(datExpr0, power = power)

TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
pdf(file=file.path(FIG_DIR, "03_gene_clustering.pdf"), width=12, height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

minModuleSize = MIN_MODULE_SIZE
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

pdf(file=file.path(FIG_DIR, "04_Dynamic_Tree.pdf"), width=8, height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# ---
# Section 7: Module Merging
# ---
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")

pdf(file=file.path(FIG_DIR, "05_Clustering_module.pdf"),width=7,height=7)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = MERGE_CUT_HEIGHT
abline(h=MEDissThres, col = "red")
dev.off()

merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

pdf(file=file.path(FIG_DIR, "06_merged_dynamic.pdf"), width = 8, height = 6)
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# ---
# Section 8: Module-Trait Correlation (Original Logic)
# ---
sameSample2=intersect(row.names(cli2), rownames(MEs))
MEs=MEs[sameSample2,]
datTraits=cli2[sameSample2,]
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

pdf(file=file.path(FIG_DIR, "07_Module_trait.pdf"), width=6.5, height=5.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# ---
# Section 9: Save Module Genes
# ---
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = file.path(TABLE_DIR, "module_all.txt"),sep="\t",row.names=F,quote=F)

for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file = file.path(MODULE_DIR, paste0("module_",modules,".txt")),
              sep="\t",row.names=F,col.names=F,quote=F)
}

# ---
# Section 10: Trait Data Re-processing (Original Logic)
# ---
datTraits=read.table("../data/input/datTrait.TXT", header=T, sep="\t", check.names=F,row.names = 1) # Assumed path
score <- read.csv("../data/input/score.csv", header = TRUE, check.names = FALSE, row.names = 1) # Assumed path
last_four_columns=score
colnames(last_four_columns) <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
datTraits <- merge(datTraits, last_four_columns, by = "row.names", all = FALSE)
rownames(datTraits) <- datTraits$Row.names
datTraits$Row.names <- NULL 

normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
datTraits_normalized <- as.data.frame(lapply(datTraits, function(col) {
  if (is.numeric(col)) {
    return(normalize(col))
  } else {
    return(col)
  }
}))
View(datTraits_normalized)
datTraits=datTraits_normalized

head(datTraits)
rm(datTraits)
datTraits1 <- datTraits[, "mRNAsi", drop = FALSE]
rownames(datTraits) <- rownames(datTraits)
datTraits=datTraits1

# ---
# Section 11: Visualize Target Module (GS vs MM)
# ---
module = "blue"
Selectedclinical = "mRNAsi"
Selectedclinical2 = "mRNAsi"

Selectedclinical = as.data.frame(datTraits[,Selectedclinical]);
names(Selectedclinical) = "Selectedclinical";
modNames = substring(names(MEs), 3)
datExpr1 = datExpr0[rownames(MEs),]
geneModuleMembership = as.data.frame(cor(datExpr1, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr1, Selectedclinical, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Selectedclinical), sep="");
names(GSPvalue) = paste("p.GS.", names(SelectedCSV), sep="");

column = match(module, modNames)
moduleGenes = moduleColors==module
outPdf=file.path(FIG_DIR, paste(Selectedclinical2,"_", module, ".pdf",sep=""))

pdf(file=outPdf,width=7,height=7)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for ", Selectedclinical2,sep=""),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# ---
# Section 12: Save Hub Genes
# ---
datMM=cbind(geneModuleMembership[,paste("MM", module, sep="")], geneTraitSignificance)
colnames(datMM)[1] = paste("MM", module, sep="")
geneSigFilter=0.9
moduleSigFilter=0.8
datMM=datMM[abs(datMM[,ncol(datMM)])>geneSigFilter,]
datMM=datMM[abs(datMM[,1])>moduleSigFilter,]
write.table(row.names(datMM), file = file.path(TABLE_DIR, paste0("hubGenes",module,".txt")),
            sep="\t",row.names=F,col.names=F,quote=F)