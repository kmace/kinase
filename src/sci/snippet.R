
library(DESeq2)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(viridis)
library(tidyr)
library(genefilter)
library(ggrepel)
library(plotly)

source('0.1-data_loading_functions.R')
source('0.2-data_manipulation_functions.R')
source('0.3-plotting_functions.R')

t2g = load_transcripts_to_genes()
sample_meta = load_sample_meta_data()
seq_meta = load_seq_meta_data()
meta = right_join(sample_meta, seq_meta, by = c("Sample_Name" = "sample"))
raw_counts = load_count_matrix(meta,'tophat')

dds_all = DESeqDataSetFromMatrix(countData = raw_counts, colData = meta, ~1)
dds_wt = dds_all[,colData(dds_all)$Strain == 'WT']

rld <- rlogTransformation(dds_wt)
rlogm <- assay(rld)

library("vsn")
library(RColorBrewer)
meanSdPlot(rlogm)
lab = colData(dds_wt)$Stress
plotHclustColors <- function(matrix,labels,hang=.1,...) {
  colnames(matrix) <- labels
  d <- dist(t(matrix))
  hc <- hclust(d)
  labelColors <- brewer.pal(nlevels(labels), "Paired")
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[which(levels(lab) == a$label)]
      attr(n, "nodePar") <- c(a$nodePar, lab.col=labCol, pch=NA)
    }
    n
  }
  clusDendro <- dendrapply(as.dendrogram(hc,hang=hang), colLab)
  plot(clusDendro,...)
}
plotHclustColors(rlogm, lab)

library(RColorBrewer)
labelColors <- brewer.pal(nlevels(lab), 'Paired')
labelColors


