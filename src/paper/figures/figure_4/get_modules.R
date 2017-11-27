library(readr)
library(igraph)
library(gtools)
library(WGCNA)
library(dplyr)

yeastnet = read_delim('../../../../input/external_datasets/YeastNet.v3.txt', delim = '\t', col_names = c('Gene_1', 'Gene_2', 'weight'))
yeastnet %>% na.omit() %>% graph.data.frame(directed = FALSE) -> g
adj = as_adjacency_matrix(g, attr='weight', type='both')
adj = as.matrix(adj)
adj_inv_logit = inv.logit(adj)
adj_inv_logit = (adj_inv_logit*2 - 1)
diag(adj_inv_logit) = 1
p_est = pickSoftThreshold.fromSimilarity(similarity = adj_inv_logit)$powerEstimate
TOM = TOMsimilarity(adjMat = adj_inv_logit^p_est)
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)


minModuleSize = 10;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);

m = data.frame(mod = dynamicMods, Gene = rownames(adj_inv_logit))
m = left_join(m, t2g)
