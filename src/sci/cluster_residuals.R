library(WGCNA)
library(dplyr)
library(tidyr)
library(caret)

load('../../input/images/paper_data.RData')

resid = t(full_std_resid_matrix)

datTraits = tibble(Experiment = rownames(resid)) %>% separate(Experiment, c('Strain_Code', 'Condition'), '_', remove = F)



resid[abs(resid)<2] = 0
datExpr0 = resid


sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#datExpr0 = resid
cutoff_height = 160
abline(h = cutoff_height, col = "red");

clust = cutreeStatic(sampleTree, cutHeight = cutoff_height, minSize = 10)

keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
datTraits = datTraits[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


sampleTree = hclust(dist(datExpr), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to check no outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


collectGarbage()

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

plot(sft_table$median.k.,
     -sign(sft_table$slope) * sft_table$SFT.R.sq,
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_table$median.k.,
     -sign(sft_table$slope) * sft_table$SFT.R.sq,
     labels=powers,
     col="red");


# should be 6, but i like 8

#sft$powerEstimate = 6

net = blockwiseModules(datExpr,
                       power = sft$powerEstimate,
                       networType = network_type,
                       # corType = 'bicor',
                       # corOptions = list(maxPOutliers = 0.05),
                       minModuleSize = 4,
                       deepSplit = 2, # 4 is very strong, 2 gets the job done. not sure about 3
                       maxBlockSize = 30000,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "kinaseTOM",
                       verbose = 3)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

table(moduleColors)
length(table(moduleColors))

# [rm(list=setdiff(ls(), c("datExpr", "sft", "network_type", "t2g", "expr", "full_std_resid_matrix", "datTraits", "resid")))

modules = tibble(Gene = colnames(datExpr),
                 module = as.character(moduleLabels),
                 color = moduleColors) %>%
                 left_join(t2g)

MEs = net$MEs;

geneTree = net$dendrograms[[1]];

nGenes = ncol(datExpr);

nSamples = nrow(datExpr);

simple_design_matrix = predict(dummyVars(~Condition + Strain, datTraits), datTraits)
#complex_design_matrix = predict(dummyVars(~Experiment, datTraits), datTraits) Waste of time, this is just a shuffled identity matrix

moduleTraitCor = cor(MEs, simple_design_matrix, use = "p");

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

library(d3heatmap)
d3heatmap(moduleTraitCor)



GOenr = GOenrichmentAnalysis(moduleLabels, yeastORFs = colnames(datExpr), organism = 'yeast')
tab = GOenr$bestPTerms[[4]]$enrichment

#all = modules %>% group_by(module) %>% nest(.key = gene_info) %>% full_join(tab %>% group_by(module) %>% nest(.key = go_info))

dir.create('WGCNA_modules')

library(readr)

modules %>%
  group_by(module) %>%
  do(write_csv(., path = paste0('WGCNA_modules/cluster_',
                                as.character(unique(.$module)),
                                '_genes_info.csv')))

write_genes = function(df, path){
write(df$name, file=path)
}

modules %>%
  group_by(module) %>%
  do(write_genes(., path = paste0('WGCNA_modules/cluster_',
                                as.character(unique(.$module)),
                                '_genes.txt')))


tab %>%
  group_by(module) %>%
  do(write_csv(., path = paste0('WGCNA_modules/cluster_',
                                as.character(unique(.$module)),
                                '_go_terms.csv')))


# After meme is run:

tfs = do.call(rbind,
lapply(dir('WGCNA_modules/tfs/'),
function(x) {
if(file.size(paste0('WGCNA_modules/tfs/', x)) > 0) {
return(data.frame(file = x, TF = read.table(paste0('WGCNA_modules/tfs/',x), stringsAsFactors = F)))
}})) %>% as_tibble %>% separate(file, c('c', 'module', 't'), sep = '_') %>% dplyr::rename(TF = V1) %>% dplyr::select(module, TF)
