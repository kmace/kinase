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
