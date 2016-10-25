# ---- Setup ----

# Master Script
library(dplyr)
library(tidyr)
library(reshape2)
library(rhdf5)
library(sleuth)
library(biomaRt)
library(DESeq2)
library(ggplot2)
library(ggbiplot)
library(gplots)
library(viridis)


select = dplyr::select
rename = dplyr::rename

source('0.1-utility_functions.R')
source('0.2-plotting_functions.R')
source('1-load_data.R')
source('2-compute_features.R')
