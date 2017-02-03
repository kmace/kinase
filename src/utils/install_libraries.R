list.of.cran.packages = c(
'devtools',
'dplyr',
'plyr',
'ggplot2',
'ggrepel',
'gplots',
'lattice',
'mclust',
'pheatmap',
'PoiClaClu',
'RColorBrewer',
'reshape2',
'tidyr',
'UpSetR',
'viridis',
'vsn',
'ggvis',
'plotly')

list.of.bioc.packages = c(
'biomaRt',
'DESeq2',
'genefilter',
'rhdf5',
'impute',
'Biobase',
'GEOquery')

new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#Poisson Distance (Witten 2011), implemented in the PoiClaClu package. This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples


source("http://bioconductor.org/biocLite.R")
biocLite()
new.packages <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

devtools::install_github("pachterlab/sleuth")
