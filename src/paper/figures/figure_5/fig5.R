library(tidyverse)
library(ComplexHeatmap)
library(circlize)
load('../../../../intermediate/images/paper_data.RData')
source('../make_obj.R')
source('../colors.R')

condition = genes %>% select(name, weights) %>% unnest() %>% filter(grepl('Condition', term)) %>% mutate(Condition = gsub(pattern =
                                                                                                                            'Condition', replacement = '', term))
strain = genes %>% select(name, weights) %>% unnest() %>% filter(grepl('Strain', term)) %>% mutate(Strain = gsub(pattern =
                                                                                                                   'Strain', replacement = '', term))

values = genes %>% select(name, data, data_augment) %>% unnest()

library(ComplexHeatmap)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm = T)
  breaks[!duplicated(breaks)]
}

# Consistnat coloring on E
col = colorRamp2(quantile_breaks(exp_mat, 11)[-c(1, 11)],
                 coolwarm_hcl)

for (target_strain in as.character(unique(values$Strain[values$Strain != 'WT']))) {
  print(target_strain)
  exp_mat = values %>% filter(Strain == target_strain) %>% select(name, Condition, Differential_Expression) %>% spread(key =
                                                                                                                         Condition, value = Differential_Expression) %>% remove_rownames() %>% column_to_rownames('name') %>% as.matrix()
  res_mat = values %>% filter(Strain == target_strain) %>% select(name, Condition, .std.resid) %>% spread(key =
                                                                                                        Condition, value = .std.resid) %>% remove_rownames() %>% column_to_rownames('name') %>% as.matrix()
  strain_weights = strain %>% filter(Strain == target_strain) %>% select(name, estimate) %>% remove_rownames() %>% column_to_rownames('name') %>% as.matrix(drop =
                                                                                                                                                              FALSE)
  strain_weights = strain_weights[rownames(exp_mat), ]
  condition_ha = HeatmapAnnotation(data.frame(Condition = colnames(exp_mat)), col = master_col)



  # ehm = Heatmap(
  #   exp_mat,
  #   #cluster_rows=F,
  #   #cluster_columns=T,
  #   col = col,
  #   show_row_names = F,
  #   show_column_names = F,
  #   use_raster = TRUE,
  #   raster_quality = 5,
  #   width = 10,
  #   #unit(4, "in"),
  #   name = 'Expression',
  #   #column_title = expression(Delta~E[ij]),
  #   row_title = expression(gene[g]),
  #   top_annotation = condition_ha
  # )
  #
  # khm = Heatmap(
  #   res_mat + strain_weights,
  #   #cluster_rows=F,
  #   #cluster_columns=T,
  #   col = col,
  #   show_row_names = F,
  #   show_column_names = F,
  #   use_raster = TRUE,
  #   raster_quality = 5,
  #   width = 10,
  #   #unit(4, "in"),
  #   name = 'Expression - K',
  #   #column_title = expression(Delta~E[ij]),
  #   row_title = expression(gene[g]),
  #   top_annotation = condition_ha
  # )
  #
  # rhm = Heatmap(
  #   res_mat,
  #   #cluster_rows=F,
  #   #cluster_columns=T,
  #   col = col,
  #   show_row_names = F,
  #   show_column_names = F,
  #   use_raster = TRUE,
  #   raster_quality = 5,
  #   width = 10,
  #   #unit(4, "in"),
  #   name = 'Residual',
  #   #column_title = expression(Delta~E[ij]),
  #   row_title = expression(gene[g]),
  #   top_annotation = condition_ha
  # )
  #
  # whm = Heatmap(
  #   strain_weights,
  #   col = col,
  #   #cluster_rows=F,
  #   #cluster_columns=T,
  #   show_row_names = F,
  #   show_column_names = F,
  #   use_raster = TRUE,
  #   raster_quality = 5,
  #   width = 2,
  #   #unit(4, "in"),
  #   name = paste('K_{', target_strain, '}'),
  #   #column_title = expression(Delta~E[ij]),
  #   row_title = expression(gene[g])#,
  #   #top_annotation = sample_ha
  # )
  # print(paste0(target_strain, '_hm.png'))
  # png(
  #   filename = paste0(target_strain, '_hm.png'),
  #   width = 6,
  #   height = 8,
  #   units = 'in',
  #   res = 300
  # )
  #
  # print(ehm  + khm + whm + rhm)
  #
  # dev.off()

#res_subset = res_mat[apply(abs(res_mat),1,max)>2, ]
res_subset = res_mat[apply(res_mat,1,function(x)diff(range(x))) > 5, ]

if(dim(res_subset)[1]>10){
library(fpc)
asw <- numeric(20)
for (k in 2:20){
  asw[[k]] <- pam(res_subset, k) $ silinfo $ avg.width
}
k.best <- which.max(asw)
} else{
  k.best=1
}

km = kmeans(res_subset,k.best)

kmdata = data.frame(Gene = names(km$cluster), cluster = km$cluster)

write_csv(kmdata, path=paste0(target_strain,'.csv'))

  rrhm = Heatmap(
    res_subset,
    #cluster_rows=F,
    #cluster_columns=T,
    col = col,
    show_row_names = F,
    show_column_names = F,
    use_raster = TRUE,
    raster_quality = 5,
    width = 10,
    #unit(4, "in"),
    name = 'Residual',
    #column_title = expression(Delta~E[ij]),
    row_title = expression(gene[g]),
    top_annotation = condition_ha,
    split=km$cluster
  )

  png(
    filename = paste0(target_strain, '_residual_small_clusters.png'),
    width = 6,
    height = 8,
    units = 'in',
    res = 300
  )

  print(rrhm)

  dev.off()
}
