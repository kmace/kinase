library(tidyverse)
library(ComplexHeatmap)
library(circlize)
load('../../../../intermediate/images/paper_data.RData')
source('../make_obj.R')
source('../colors.R')









quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=T)
  breaks[!duplicated(breaks)]
}

make_hm = function(mat, ...){
    hm = Heatmap(mat,
                 col = colorRamp2(
                        quantile_breaks(mat, 11)[-c(1, 11)],
                        coolwarm_hcl),
        show_row_names=F,
        show_column_names=F,
        ...)
    return(hm)

}

sample_ha =  HeatmapAnnotation(sample_meta,
                               col = master_col,
                               show_annotation_name = TRUE,
                               annotation_legend_param = list(Kinase = list(ncol = 1,
                                                                            title_position = "topleft",
                                                                            by_row=FALSE))
                                                                        )


e_hm = make_hm(exp_matrix,
    name = '',
    cluster_columns = F,
    top_annotation = sample_ha,
    heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
    use_raster = TRUE,
    raster_quality = 4)



pdf('heatmap.pdf', width = 6, height = 8)
draw(e_hm, heatmap_legend_side = "bottom")#, annotation_legend_side='bottom')
dev.off()

png('heatmap.png', res = 300, width = 6, height = 8)
draw(e_hm, heatmap_legend_side = "bottom")
dev.off()


library(ggrepel)
library(ggthemes)

pdf('tsne.pdf', width=5.5, height=5.5)
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
