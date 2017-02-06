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


# Example of Tunicamycin Plot
x = assay(vsd)
x = x - rowMeans(x) # Can we do this in a better way?
colnames(x) = meta$Sample_Name
t = vsd[,colData(vsd)$Stress == 'Tunicamycin']
tx = assay(t)
tm = colData(t)
rm = rowMedians(tx[,tm$Strain=='WT'])
dec = order(rm, decreasing = T)
tx = tx[dec,] # Can we also subset these so that we only keep genes that are actually on? what about DE?
rm = rm[dec]
colnames(tx) = tm$Strain_Code

sd = apply(msm,1,sd)

sigma = apply(tx,2,function(x) (x-rm)/sd)
sigma = rename_gene_names(sigma,t2g = t2g)
sigma[sigma < 1 & sigma > -1] = 0

pheatmap(sigma[1:30,]) # again, lets have this be subset, not some dumb number
pheatmap(sigma[1:300,])
pheatmap(sigma[1:1000,])


# write to GCT for morpheus
sgn = data.frame(sigma)
t = t2g[match(rownames(sgn),t2g$name),]
tmm = tm[,1:9]
writeGCT(sgn, t, tmm, 'Tunicamycin_2_sigma.gct')
