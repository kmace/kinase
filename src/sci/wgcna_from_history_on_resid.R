library(WGCNA)
library(dplyr)
library(tidyr)
library(caret)

load('input/images/paper_data.RData')

expr = t(full_std_resid_matrix)

datTraits = tibble(Experiment = rownames(expr)) %>% separate(Experiment, c('Condition', 'Strain'), '_', remove = F)



expr[abs(expr)<2] = 0
datExpr0 = expr


sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#datExpr0 = expr
cutoff_height = 160
abline(h = cutoff_height, col = "red");

clust = cutreeStatic(sampleTree, cutHeight = cutoff_height, minSize = 10)

keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datTraits = datTraits[keepSamples, ]

sampleTree = hclust(dist(datExpr), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to check no outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


collectGarbage();

enableWGCNAThreads()




network_type = 'signed'


powers = c(1:40)

sft = pickSoftThreshold(datExpr,
                        powerVector = powers,
                        verbose = 5,
                        # corFnc = bicor,
                        # corOptions = list(maxPOutliers = 0.05),
                        networkType = network_type)

sft_table = sft$fitIndices

plot(sft_table$Power,
     -sign(sft_table$slope) * sft_table$SFT.R.sq,
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_table$Power,
     -sign(sft_table$slope) * sft_table$SFT.R.sq,
     labels=powers,
     col="red");


# should be 10, but ive been using 7



net = blockwiseModules(datExpr,
                       power = sft$powerEstimate,
                       TOMType = network_type,
                       # corType = 'bicor',
                       # corOptions = list(maxPOutliers = 0.05),
                       minModuleSize = 8,
                       reassignThreshold = 0,
                       maxBlockSize = 50000,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "kinaseTOM",
                       verbose = 3)

mergedColors = labels2colors(net$colors)

moduleLabels = net$colors


MEs = net$MEs;

geneTree = net$dendrograms[[1]];

nGenes = ncol(datExpr);

nSamples = nrow(datExpr);

simple_design_matrix = predict(dummyVars(~Condition + Strain, datTraits), datTraits)
complex_design_matrix = predict(dummyVars(~Experiment, datTraits), datTraits)

moduleTraitCor = cor(MEs, complex_design_matrix, use = "p");

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

library(d3heatmap)
d3heatmap(moduleTraitCor)



GOenr = GOenrichmentAnalysis(moduleLabels, yeastORFs = colnames(datExpr), organism = 'yeast')
tab = GOenr$bestPTerms[[4]]$enrichment

modules = tibble(Gene = colnames(datExpr), module = as.character(moduleLabels)) %>% left_join(t2g) %>% group_by(module) %>% nest(.key = gene_info) %>% full_join(tab %>% group_by(module) %>% nest(.key = go_info))

print_module_genes = function(df){

}
