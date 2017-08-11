# parameters %>% filter(name == 'HSP12') %>% View()
#
# parameters %>% dplyr::filter(abs(estimate)>1 & p.value < 0.01) %>% group_by(term, model_type) %>% summarise(num_genes = n())
#
# parameters %>% dplyr::filter(abs(estimate)>1 & p.value < 0.01) %>% spread(Gene,term,estimate)
#
# good_parameters = parameters %>% dplyr::filter(abs(estimate)>1 & p.value < 0.01)
# lists = sapply(terms, function(x) select(filter(good_parameters,term==x),Gene))
#
# estimates = parameters %>% select(Gene, term, estimate) %>% spread(term, estimate) %>% ungroup()
# rownames(estimates) = estimates$Gene
#
# estimates_matrix = as.matrix(select(estimates, -Gene))
# col_mean = apply(estimates_matrix,2,mean)
# #m = apply(parameter_matrix,1,function(x) x - col_m)
# #m = t(m)
# #heatmap3(m[m[,10]>1,-10], labRow = NA, scale='none')
# heatmap3(parameter_matrix, labRow = NA, scale='none')
#
#
# library(plotly)
# library(d3heatmap)
#
# shrink_large = function(data, max, min = -max){
#   data[data>max] = max
#   data[data<min] = min
#   return(data)
# }
#
# d3heatmap(a[str_length(rownames(a)) < 7,grep('Cond', colnames(a))], Rowv = NA, Colv=NA)
#
# heatmaply(shrink_large(a[str_length(rownames(a)) < 7,
#                          grep('Cond', colnames(a))], 2),
#           Rowv = NA,
#           Colv=NA,
#           scale_fill_gradient_fun= scale_fill_gradient2())