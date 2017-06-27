library(heatmap3)
library(gputools)

rm(list=ls())
load('../../input/images/normalized_data.RData')
sds = apply(vlog,1,sd)
colnames(vlog) = meta$Experiment
vlog = vlog[sds>.2,]
gene_cor = cor(t(vlog))
gene_cor_dist = as.dist(abs(gene_cor - 1))
gene_cor_dist_hclust = hclust(gene_cor_dist)
gene_cor_ordered = gene_cor[gene_cor_dist_hclust$order, gene_cor_dist_hclust$order]


pdf('gene_heatmap.pdf')
heatmap3(gene_cor_ordered,
         Colv = NA,
         Rowv = NA,
         labCol = NA,
         labRow = NA,
         symm=TRUE,
         scale = 'none')
dev.off()


ann = meta[,c("Condition"), drop=F]
ann_colors = cbind(b[ann$Condition])
colnames(ann_colors) = colnames(ann)
pdf('gene_heatmap_ordered.pdf')
heatmap3(vlog[o,],
Colv =NA,
Rowv = NA,
labCol = NA,
labRow = NA,
scale = 'none',
ColSideColors = ann_colors,
useRaster = F)
dev.off()

pdf('all_heatmap.pdf')
heatmap3(vlog_ordered,
         Colv = NA,
         Rowv = NA,
         labCol = NA,
         labRow = NA,
         scale = 'row')
dev.off()

# gene_cor_squashed = gene_cor
# gene_cor_squashed[abs(gene_cor)<.6] = 0
# 
# gene_cor_dist = as.dist(abs(gene_cor_squashed - 1))
# gene_cor_dist_hclust = hclust(gene_cor_dist)
# gene_cor_ordered = gene_cor[gene_cor_dist_hclust$order, gene_cor_dist_hclust$order]
# 
# library(heatmap3)
# pdf('gene_heatmap.pdf')
# heatmap3(gene_cor_ordered,
#          Colv = NA,
#          Rowv = NA,
#          labCol = NA,
#          labRow = NA,
#          scale = 'none')
# dev.off()


myheatmap = function(x) {
  heatmap3(x,
Colv = NA,
Rowv = NA,
labCol = NA,
labRow = NA,
scale = 'none')
}
