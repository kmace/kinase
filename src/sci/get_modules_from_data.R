library(readr)
library(igraph)
library(gtools)
library(WGCNA)
library(dplyr)
library(tidyr)
library(caret)
library(tidyverse)

get_module_from_graph = function(g){

    adj = as_adjacency_matrix(g, attr='weight', type='both')
    adj = as.matrix(adj)
    adj_inv_logit = inv.logit(adj)
    adj_inv_logit = (adj_inv_logit * 2 - 1)
    diag(adj_inv_logit) = 1
    p_est = pickSoftThreshold.fromSimilarity(similarity = adj_inv_logit)$powerEstimate
    TOM = TOMsimilarity(adjMat = adj_inv_logit^p_est)
    dissTOM = 1-TOM
    # Call the hierarchical clustering function
    geneTree = hclust(as.dist(dissTOM), method = "average");
    # Plot the resulting clustering tree (dendrogram)
    minModuleSize = 15;
    # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = geneTree,
                                distM = dissTOM,
                                deepSplit = 2,
                                pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize);

    modules = data.frame(module = dynamicMods,
                         Gene = rownames(adj_inv_logit)) %>% left_join(t2g)
    return(modules)
}

filter_samples = function(exp){
    sampleTree = hclust(dist(exp), method = "average");
    #sizeGrWindow(12,9)
    #par(cex = 0.6);
    #par(mar = c(0,4,2,0))
    #plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
    #     cex.axis = 1.5, cex.main = 2)

    #datExpr0 = resid
    cutoff_height = 160
    #abline(h = cutoff_height, col = "red"); # Check for the right cutoff

    clust = cutreeStatic(sampleTree, cutHeight = cutoff_height, minSize = 15)

    keepSamples = (clust==1)

    return(keepSamples)
}

get_module_from_exp = function(datExpr){
    network_type = 'signed'
    powers = c(1:40)
    sft = pickSoftThreshold(datExpr,
                            powerVector = powers,
                            verbose = 5,
                            # corFnc = bicor,
                            # corOptions = list(maxPOutliers = 0.05),
                            networkType = network_type)

    sft_table = sft$fitIndices
    print(names(sft))
    net = blockwiseModules(datExpr,
                           power = sft$powerEstimate,
                           networType = network_type,
                           # corType = 'bicor',
                           # corOptions = list(maxPOutliers = 0.05),
                           minModuleSize = 15,
                           deepSplit = 2, # 4 is very strong, 2 gets the job done. not sure about 3
                           #maxBlockSize = 30000,
                           numericLabels = TRUE,
                           pamRespectsDendro = FALSE,
                           #saveTOMs = TRUE,
                           #saveTOMFileBase = "kinaseTOM",
                           verbose = 3)

    moduleLabels = net$colors
    moduleColors = labels2colors(moduleLabels)


    modules = tibble(Gene = colnames(datExpr),
                     module = as.character(moduleLabels),
                     color = moduleColors) %>%
                     left_join(t2g)

return(modules)

    # adj = as_adjacency_matrix(g, attr='weight', type='both')
    # adj = as.matrix(adj)
    # adj_inv_logit = inv.logit(adj)
    # adj_inv_logit = (adj_inv_logit * 2 - 1)
    # diag(adj_inv_logit) = 1
    # p_est = pickSoftThreshold.fromSimilarity(similarity = adj_inv_logit)$powerEstimate
    # TOM = TOMsimilarity(adjMat = adj_inv_logit^p_est)
    # dissTOM = 1-TOM
    # # Call the hierarchical clustering function
    # geneTree = hclust(as.dist(dissTOM), method = "average");
    # # Plot the resulting clustering tree (dendrogram)
    # minModuleSize = 15;
    # # Module identification using dynamic tree cut:
    # dynamicMods = cutreeDynamic(dendro = geneTree,
    #                             distM = dissTOM,
    #                             deepSplit = 2,
    #                             pamRespectsDendro = FALSE,
    #                             minClusterSize = minModuleSize);
    #
    # modules = data.frame(module = dynamicMods,
    #                      Gene = rownames(adj_inv_logit)) %>% left_join(t2g)
    # return(modules)
}

load('intermediate/images/paper_data.RData')

enableWGCNAThreads()


yeastnet = read_delim('input/external_datasets/module_definitions/YeastNet.v3.txt', delim = '\t', col_names = c('Gene_1', 'Gene_2', 'weight'))
yeastnet %>% na.omit() %>% graph.data.frame(directed = FALSE) -> g
yeastnet_modules = get_module_from_graph(g)



############ yeastnet above, my data below


exp = genes %>%
    select(data, name) %>%
    unnest() %>%
    select(name, Sample_Name, Expression) %>%
    spread(key=name, value=Expression) %>%
    remove_rownames() %>%
    column_to_rownames('Sample_Name') %>%
    as.matrix()

# what about using the corrolation of the outlier probabilites as the adjacency matrix? WGCNA does encorage use of datExpr though...

resid = genes %>%
    select(data, data_augment, name) %>%
    unnest() %>%
    select(name, Sample_Name, .resid) %>%
    spread(key=name, value=.resid) %>%
    remove_rownames() %>%
    column_to_rownames('Sample_Name') %>%
    as.matrix()

#resid[abs(resid)<2] = 0
sample_keep = filter_samples(resid)
resid = resid[sample_keep,]
resid_modules = get_module_from_exp(resid)[,1:3]

sample_keep = filter_samples(exp)
exp = exp[sample_keep,]
exp_modules = get_module_from_exp(exp)[,1:3]


colnames(resid_modules)[1] = 'name'
colnames(exp_modules)[1] = 'name'

load(file = "intermediate/images/externally_defined_modules.RData")
save(list = ls()[grep('module', ls())], file = "intermediate/images/externally_defined_modules.RData")

