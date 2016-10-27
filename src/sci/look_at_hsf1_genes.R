library(DESeq2)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(viridis)
library(tidyr)
library(genefilter)
library(vsn)
library(RColorBrewer)


source('../utils/0.1-data_loading_functions.R')
source('../utils/0.2-data_manipulation_functions.R')
source('../utils/0.3-plotting_functions.R')

load('../../input/images/normalized_data.RData')
#hsf1_genes = read.table('../../input/reference/HSF1_targets.txt', header=F, stringsAsFactors=F)[,1]
HS_genes = read.table('../../input/reference/Heatshock_Interest_targets.txt', header=F, stringsAsFactors=F)[,1]
# Normalize
control_drug = which(
    colData(rld_all)$Stress == 'None' &
    colData(rld_all)$Drug == 'Cocktail' &
    colData(rld_all)$Media == 'YPD' &
    colData(rld_all)$Strain == 'WT')

control_mean = apply((assay(rld_all[,control_drug])),1,mean)
norm = apply(assay(rld_all),2, function(x) x - control_mean)

# Get Genes
HS_idx = which(t2g$name %in% HS_genes)
HS_target_id = t2g$target_id[HS_idx]

# Get Samples
heatshock_drug_idx = which(colData(rld_all)$Stress == 'Heatshock' & colData(rld_all)$Drug == 'Cocktail')
salt_drug_idx = which(colData(rld_all)$Stress == 'Salt' & colData(rld_all)$Drug == 'Cocktail')

subset = norm[HS_target_id, heatshock_drug_idx]

# Rename cols and rows
colnames(subset) = colData(rld_all)$Strain[heatshock_drug_idx]
rownames(subset) = t2g$name[hsf1_idx]
