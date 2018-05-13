library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(outliers)
load('../../../../intermediate/images/paper_data.RData')
source('../make_obj.R')
source('../colors.R')
load('../../../../intermediate/images/externally_defined_modules.RData')
condition = genes %>% select(name, weights) %>% unnest() %>% filter(grepl('Condition', term)) %>% mutate(Condition = gsub(pattern =
                                                                                                                            'Condition', replacement = '', term))
strain = genes %>% select(name, weights) %>% unnest() %>% filter(grepl('Strain', term)) %>% mutate(Strain = gsub(pattern =
                                                                                                                   'Strain', replacement = '', term))

values = genes %>% select(name, data, data_augment) %>% unnest() %>% left_join(t2g)


library(ComplexHeatmap)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm = T)
  breaks[!duplicated(breaks)]
}

# Consistnat coloring on E
col = colorRamp2(quantile_breaks(exp_matrix, 11)[-c(1, 11)],
                 coolwarm_hcl)

for (target_cond in as.character(unique(values$Condition[values$Condition != 'WT']))) {

  #target_cond = 'Tunicamycin'
  mod_values  = values %>% left_join(modules) %>% group_by(Strain_Code, Condition, module) %>% summarize_if(is.numeric,mean, na.rm=T) %>% select(module, .std.resid) %>% ungroup() %>% left_join(modules %>% group_by(module) %>% summarize(size = n())) %>% mutate(module_size_corrected_resid = .std.resid * size)
  mod_values %>% select(Condition, Strain_Code, module, .std.resid) %>% na.omit() %>% filter(Condition == target_cond) %>% spread(key=Strain_Code, value = .std.resid) %>% select_if(is.numeric) %>% remove_rownames() %>% column_to_rownames('module') %>% as.matrix() -> h
  #heatmaply(h, colors =  PiYG)
  h = h - apply(h[,grepl('WT',colnames(h))],1,mean)

  #heatmaply(h,colors = PiYG)

  rrhm = Heatmap(
    h,
    show_row_names = F,
    show_column_names = F,
    use_raster = TRUE,
    raster_quality = 5)

  pdf(paste0(target_cond,'_heatmap.pdf'))
  print(rrhm)
  dev.off()


  con = abs(h) > 2
  g <- graph.incidence(t(con), directed = T, mode = 'out')
  gsub = igraph::delete.vertices(igraph::simplify(g), igraph::degree(g)==0)

  pdf(paste0(target_cond,'_graph.pdf'))
  plot(gsub, vertex.color = c('pink','white','yellow','red')[2+sign(igraph::degree(gsub,mode='out') - igraph::degree(gsub,mode='in'))],vertex.label.cex = .5, vertex.size = 5, edge.arrow.size=.3, edge.color = 'black', layout = layout_nicely(gsub, niter = 100000)*2, main=target_cond)
  dev.off()
  #   print(target_strain)
  #   exp_mat = values %>%
  #     filter(Condition == target_strain) %>%
  #     select(name, Strain, Differential_Expression) %>%
  #     spread(key = Strain, value = Differential_Expression) %>%
  #     remove_rownames() %>% column_to_rownames('name') %>%
  #     as.matrix()
  #
  #   res_mat = values %>%
  #     filter(Condition == target_condition) %>%
  #     select(name, Strain_Code, .std.resid) %>%
  #     spread(key = Strain_Code, value = .std.resid) %>%
  #     remove_rownames() %>% column_to_rownames('name') %>%
  #     as.matrix()
  #
  #   strain_weights = strain %>% filter(Strain == target_strain) %>% select(name, estimate) %>% remove_rownames() %>% column_to_rownames('name') %>% as.matrix(drop =
  #                                                                                                                                                               FALSE)
  #   strain_weights = strain_weights[rownames(exp_mat), ]
  #   condition_ha = HeatmapAnnotation(data.frame(Strain = colnames(res_mat)), col = master_col)
  #
  #
  #
  #   # ehm = Heatmap(
  #   #   exp_mat,
  #   #   #cluster_rows=F,
  #   #   #cluster_columns=T,
  #   #   col = col,
  #   #   show_row_names = F,
  #   #   show_column_names = F,
  #   #   use_raster = TRUE,
  #   #   raster_quality = 5,
  #   #   width = 10,
  #   #   #unit(4, "in"),
  #   #   name = 'Expression',
  #   #   #column_title = expression(Delta~E[ij]),
  #   #   row_title = expression(gene[g]),
  #   #   top_annotation = condition_ha
  #   # )
  #   #
  #   # khm = Heatmap(
  #   #   res_mat + strain_weights,
  #   #   #cluster_rows=F,
  #   #   #cluster_columns=T,
  #   #   col = col,
  #   #   show_row_names = F,
  #   #   show_column_names = F,
  #   #   use_raster = TRUE,
  #   #   raster_quality = 5,
  #   #   width = 10,
  #   #   #unit(4, "in"),
  #   #   name = 'Expression - K',
  #   #   #column_title = expression(Delta~E[ij]),
  #   #   row_title = expression(gene[g]),
  #   #   top_annotation = condition_ha
  #   # )
  #   #
  #   # rhm = Heatmap(
  #   #   res_mat,
  #   #   #cluster_rows=F,
  #   #   #cluster_columns=T,
  #   #   col = col,
  #   #   show_row_names = F,
  #   #   show_column_names = F,
  #   #   use_raster = TRUE,
  #   #   raster_quality = 5,
  #   #   width = 10,
  #   #   #unit(4, "in"),
  #   #   name = 'Residual',
  #   #   #column_title = expression(Delta~E[ij]),
  #   #   row_title = expression(gene[g]),
  #   #   top_annotation = condition_ha
  #   # )
  #   #
  #   # whm = Heatmap(
  #   #   strain_weights,
  #   #   col = col,
  #   #   #cluster_rows=F,
  #   #   #cluster_columns=T,
  #   #   show_row_names = F,
  #   #   show_column_names = F,
  #   #   use_raster = TRUE,
  #   #   raster_quality = 5,
  #   #   width = 2,
  #   #   #unit(4, "in"),
  #   #   name = paste('K_{', target_strain, '}'),
  #   #   #column_title = expression(Delta~E[ij]),
  #   #   row_title = expression(gene[g])#,
  #   #   #top_annotation = sample_ha
  #   # )
  #   # print(paste0(target_strain, '_hm.png'))
  #   # png(
  #   #   filename = paste0(target_strain, '_hm.png'),
  #   #   width = 6,
  #   #   height = 8,
  #   #   units = 'in',
  #   #   res = 300
  #   # )
  #   #
  #   # print(ehm  + khm + whm + rhm)
  #   #
  #   # dev.off()
  #
  # #res_subset = res_mat[apply(abs(res_mat),1,max)>2, ]
  # res_subset = res_mat[apply(res_mat,1,function(x)diff(range(x))) > 5, ]
  #
  # if(dim(res_subset)[1]>10){
  # library(fpc)
  # asw <- numeric(20)
  # for (k in 2:20){
  #   asw[[k]] <- pam(res_subset, k) $ silinfo $ avg.width
  # }
  # k.best <- which.max(asw)
  # } else{
  #   k.best=1
  # }
  #
  # km = kmeans(res_subset,k.best)
  #
  # kmdata = data.frame(Gene = names(km$cluster), cluster = km$cluster)
  #
  # write_csv(kmdata, path=paste0(target_strain,'.csv'))
  #
  #   rrhm = Heatmap(
  #     res_subset,
  #     #cluster_rows=F,
  #     #cluster_columns=T,
  #     col = col,
  #     show_row_names = F,
  #     show_column_names = F,
  #     use_raster = TRUE,
  #     raster_quality = 5,
  #     width = 10,
  #     #unit(4, "in"),
  #     name = 'Residual',
  #     #column_title = expression(Delta~E[ij]),
  #     row_title = expression(gene[g]),
  #     top_annotation = condition_ha,
  #     split=km$cluster
  #   )
  #
  #   png(
  #     filename = paste0(target_strain, '_residual_small_clusters.png'),
  #     width = 6,
  #     height = 8,
  #     units = 'in',
  #     res = 300
  #   )
  #
  #   print(rrhm)

  #  dev.off()
}



    #biggest_dist = lapply(as.numeric(names(table(table(m$mod)))), function(n) replicate(10000, max(abs(apply(std_resid_matrix[sample(1:dim(std_resid_matrix)[1], n, FALSE),],2,mean))))); names(biggest_dist) = names(table(table(m$mod)))
#go_slim = read_delim('../../../../input/external_datasets/go_slim_mapping.tab', '\t', col_names = F)
#colnames(go_slim) = c('Gene', 'name', 'sdg', 'go_class', 'term', 'go_num', 'geen_type')
#go_values  = values %>% left_join(go_slim) %>% group_by(Strain_Code, Condition, term) %>% summarize_if(is.numeric,mean, na.rm=T) %>% select(term, .std.resid) %>% ungroup() %>% left_join(go_slim %>% group_by(term) %>% summarize(size = n()))
#go_values %>% select(Condition, Strain_Code, term, .std.resid) %>% na.omit() %>% filter(Condition == target_cond) %>% spread(key=Strain_Code, value = .std.resid) %>% remove_rownames() %>% column_to_rownames('term') %>% select_if(is.numeric) %>% as.matrix() -> h
#heatmaply::heatmaply(h^7)

#library(data.table)
#reg = load_binary_yeastract_regulation_table('../../../../input/external_datasets/yeastract/RegulationTwoColumnTable_Documented_both_evidence_stress.tsv')
#reg = as.data.frame(reg)
#colnames(reg) = c('TF', 'name')
#vals = values %>% select(Condition, name, Strain_Code, .std.resid, .fitted, Expression)
#tf_val = full_join(vals, reg)
#tf_res = tf_val %>% group_by(Condition, Strain_Code, TF) %>% summarize_if(is.numeric,mean) %>% ungroup() %>% arrange(desc(abs(.std.resid)))
#hist(tf_res$.std.resid,100)
#tf_res %>% filter(abs(.std.resid)>2)
library(ggnetwork)
library(network)
library(sna)
#network(tf_res %>% filter(abs(.std.resid)>2 & Condition == "YPD"))
ns = names(modules)

for (i in seq_along(modules)){
  mod = modules[[i]]
  mod_res = values %>%
    full_join(mod) %>%
    select(Condition, name, Strain_Code, .std.resid, .fitted, Expression, module) %>%
    group_by(Condition, Strain_Code, module) %>%
    summarize_if(is.numeric,mean) %>%
    ungroup() %>%
    arrange(desc(abs(.std.resid)))


  pdf(paste0(ns[[i]], '_networks.pdf'))
  for (target_cond in as.character(unique(values$Condition))) {
  #tf_res %>% filter(Condition == target_cond) %>% select(-Condition,-.fitted,-Expression) %>% spread(key=Strain_Code, value = .std.resid) %>% na.omit() %>% remove_rownames() %>% column_to_rownames('TF') %>% as.matrix() %>% apply(1,function(x) x - mean(x[grepl('WT',names(x))])) %>% apply(2,function(x) scores(x, prob=.99)) %>% igraph::graph.incidence(directed = TRUE, mode = 'out') -> g
  #tf_exp %>% filter(Condition == target_cond) %>% select(-Condition) %>% spread(key=Strain_Code, value = Exp) %>% na.omit() %>% remove_rownames() %>% column_to_rownames('TF') %>% as.matrix() %>% t() %>% apply(2,function(x) scores(x, prob=.98)) %>% igraph::graph.incidence(directed = TRUE, mode = 'out') -> g
  #tf_exp %>% filter(Condition == target_cond) %>% select(-Condition) %>% spread(key=Strain_Code, value = Exp) %>% na.omit() %>% remove_rownames() %>% column_to_rownames('TF') %>% as.matrix() %>% apply(1,function(x) x - mean(x[grepl('WT',names(x))])) %>% apply(2,function(x) scores(x, prob=.99)) %>% igraph::graph.incidence(directed = TRUE, mode = 'out') -> g
  # gsub = igraph::delete.vertices(igraph::simplify(g), igraph::degree(g)==0)
  #segal_res %>% filter(Condition == target_cond) %>% select(-Condition,-.std.resid) %>% spread(key=Strain_Code, value = Expression) %>% na.omit() %>% remove_rownames() %>% column_to_rownames('Module Description') %>% select(-module) %>% select_if(is.numeric) %>% as.matrix() %>% apply(1,function(x) x - mean(x[grepl('WT',names(x))])) %>% apply(2,function(x) scores(x, prob=.99)) %>% igraph::graph.incidence(directed = TRUE, mode = 'out') -> g
  mod_res %>%
      filter(Condition == target_cond) %>%
      select(Strain_Code,module,.std.resid) %>%
      spread(key=Strain_Code, value = .std.resid) %>%
      na.omit() %>%
      remove_rownames() %>%
      column_to_rownames('module') %>%
      select_if(is.numeric) %>%
      as.matrix() %>%
      apply(1,function(x) x - mean(x[grepl('WT',names(x))])) %>%
      apply(2,function(x) scores(x, prob=.99)) %>%
      igraph::graph.incidence(directed = TRUE, mode = 'out') -> g
  gsub = igraph::delete.vertices(igraph::simplify(g), igraph::degree(g)==0)
  if(length(E(gsub))>0) plot(gsub, vertex.color = c('pink','white','yellow','red')[2+sign(igraph::degree(gsub,mode='out') - igraph::degree(gsub,mode='in'))],vertex.label.cex = .5, vertex.size = 5, edge.arrow.size=.3, edge.color = 'black', layout = layout_nicely(gsub, niter = 100000)*2, main=target_cond)
}
dev.off()
}





for (target_cond in as.character(unique(values$Condition))) {
  library(graph)
  library(RCy3)

  dat = mod_res %>%
    filter(Condition == target_cond) %>%
    select(-Condition) %>% na.omit()

  V = na.omit(unique(c(dat$module,dat$Strain_Code)))

  edLM = dat %>%
    filter(abs(.std.resid)>.5)  %>%
    split(.$Strain_Code) %>%
    map(~ .x %>% select(-Strain_Code) %>%
          rename(edges = module, weights = Expression) %>%
          mutate(edges = match(edges, V)))
  edL <- vector("list", length=length(V))
  names(edL) = V

  for (i in 1:length(edL)){
    edL[[i]] = list(edges=character(),weights=numeric())
  }
  for (i in 1:length(edLM)){
    edL[[names(edLM)[i]]] = list(edges=edLM[[i]]$edges,weights=edLM[[i]]$weights)
  }

  g = graphNEL(nodes=V, edgeL=edL, edgemode = 'directed')
  g <- initNodeAttribute (graph=g,
                          attribute.name='moleculeType',
                          attribute.type='char',
                          default.value='undefined')

  g <- initEdgeAttribute (graph=g,
                          attribute.name='weight',
                          attribute.type='numeric',
                          default.value='0')

  for(kin in unique(dat$Strain_Code)){
    nodeData (g, kin, 'moleculeType') <- 'Kinase'
  }

  for(tf in unique(dat$module)){
    nodeData (g, tf, 'moleculeType') <- 'Module'
  }

  cw <- CytoscapeWindow (target_cond, graph=g, overwrite=TRUE)
  displayGraph (cw)
  setVisualStyle('Curved')
  ?setLayoutProperties
}
