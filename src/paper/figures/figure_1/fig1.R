library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(ggthemes)

# Objects Required:
# sample_meta
# master_col
# exp_matrix
# sample_tsne
# condition_color_scale


load('../../../../intermediate/images/paper_data.RData')
source('../make_obj.R')
source('../colors.R')


#Set output location:
output_path = '../../../../output/Images/figure_1'
dir.create(output_path)
# Heatmap

#quantile_breaks <- function(xs, n = 10) {
#  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=T)
#  breaks[!duplicated(breaks)]
#}

# Column Annotation
column_annotation =  HeatmapAnnotation(sample_meta,
                               col = master_col,
                               show_annotation_name = TRUE,
                               annotation_legend_param = list(Kinase = list(ncol = 2,
                                                                            title_position = "topleft",
                                                                            by_row=FALSE))
                                                                        )
# Cluster Rows
split = cutree(hclust(dist(exp_matrix),
                      method = 'ward.D2'),
               k=10)

# Heatmap
hmap = Heatmap(exp_matrix,
               name = '',
               cluster_columns = F,
               top_annotation = column_annotation,
               heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
               use_raster = FALSE,
               #col = colorRamp2(quantile_breaks(exp_matrix, 11), divergent_colors),
               col = colorRamp2(-5:5, divergent_colors),
               split = split,
               show_row_names=F,
               show_column_names=F)





pdf(file.path(output_path, 'heatmap.pdf'), width = 16, height = 9)
draw(hmap, heatmap_legend_side = "bottom")#, annotation_legend_side='bottom')
dev.off()

#png(file.path(output_path, 'heatmap.png'), res = 300, width = 16, height = 9)
#draw(hmap, heatmap_legend_side = "bottom")
#dev.off()



pdf(file.path(output_path, 'tsne.pdf'), width=5.5, height=5.5)
ggplot(sample_tsne,
       aes(x=TSNE1,
           y=TSNE2,
           color = Condition,
           label = Strain_Code)) +
  geom_point(size=4, alpha=.6) +
  geom_text_repel(size = 1.5,
                  color = 'black',
                  point.padding = NA,
                  box.padding = unit(0.01, "lines")) +
  theme_Publication() +
  theme(legend.position="none") +
  condition_color_scale
dev.off()


# theme_fivethirtyeight(base_family='') +
# theme(plot.background = element_blank(),
#       legend.background = element_blank(),
#       legend.key = element_blank(),
#       panel.background = element_blank()) +
