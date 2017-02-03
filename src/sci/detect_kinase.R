get_extream_kinase = function(data, strain_code, num_top, get_highest_kinase=TRUE, expression_high_to_low = TRUE) {
  data_mean = rowMeans(data)
  data = data[order(data_mean, decreasing = expression_high_to_low),]
  id = apply(data,1,function(x) strain_code[order(x,decreasing = get_highest_kinase)][1:num_top])
  dim(id) <- NULL
  rank_df = data.frame(id,count=ave(id==id, id, FUN=cumsum))
  rank_order = rep(1:dim(data)[1], each=num_top)
  rank_df$rank_order = rank_order
  return(rank_df)
}


make_plot = function(data, meta, title){
  data = data[order(rowMeans(data)), ]
   top = get_extream_kinase(data, 
                            meta$Strain_Code, 
                            num_top =3, 
                            get_highest_kinase=FALSE, # Lowest
                            expression_high_to_low = TRUE)
   plot = ggplot(top, aes(x=rank_order, y=count, group=id, color=id)) + geom_line() + labs(title = title) 
   ggplotly(plot) 
}













library(ggplot2)
vsd = vsd[,colData(vsd)$Drug == 'Cocktail']
meta = colData(vsd)
meta = as.data.frame(meta)
meta = meta %>% mutate(Condition = if_else(Stress=='None',
                                           as.character(Media),
                                           as.character(Stress)
))

meta$Condition = as.factor(meta$Condition)
x = assay(vsd)
x = x - rowMeans(x)
colnames(x) = meta$Sample_Name
data = x

for(c in unique(meta$Condition)){
  temp_idx = which(meta$Condition == c)
  temp_data = data[,temp_idx]
  temp_meta = meta[temp_idx,]
  make_plot(temp_data, temp_meta, title=c)
}