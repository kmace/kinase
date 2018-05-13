g = graphNEL(nodes=V, edgeL=edL, edgemode = 'directed')
g <- initNodeAttribute (graph=g,  attribute.name='moleculeType',
attribute.type='char',
default.value='undefined')
g <- initEdgeAttribute (graph=g,  attribute.name='weight',
attribute.type='numeric',
default.value='0')
for(kin in unique(dat$Strain_Code)){
nodeData (g, kin, 'moleculeType') <- 'Kinase'
}
for(tf in unique(dat$TF)){
nodeData (g, tf, 'moleculeType') <- 'TF'
}
cw <- CytoscapeWindow (target_cond, graph=g, overwrite=TRUE)
displayGraph (cw)
layoutNetwork(cw, 'circular')
V = unique(c(dat$TF,dat$Strain_Code))
edL = dat %>% split(.$Strain_Code) %>% map(~ .x %>% select(-Strain_Code) %>% rename(edges = TF, weights =Exp))
dat = tf_exp %>% filter(Condition == target_cond) %>% select(-Condition)


dat = mod_res %>% filter(Condition == target_cond) %>% select(-Condition) %>% na.omit()
V = na.omit(unique(c(dat$module,dat$Strain_Code)))
edLM = dat %>% filter(Expression>1)  %>% split(.$Strain_Code) %>% map(~ .x %>% select(-Strain_Code) %>% rename(edges = module, weights = Expression) %>% mutate(edges = match(edges, V)))
edL <- vector("list", length=length(V))
names(edL) = V
for (i in 1:length(edL)){edL[[i]] = list(edges=character(),weights=numeric())}
for (i in 1:length(edLM)){edL[[names(edLM)[i]]] = list(edges=edLM[[i]]$edges,weights=edLM[[i]]$weights)}
g = graphNEL(nodes=V, edgeL=edL, edgemode = 'directed')
g <- initNodeAttribute (graph=g,  attribute.name='moleculeType',
g <- initEdgeAttribute (graph=g,  attribute.name='weight',
attribute.type='numeric',
default.value='0')
for(kin in unique(dat$Strain_Code)){
nodeData (g, kin, 'moleculeType') <- 'Kinase'
}
for(tf in unique(dat$TF)){
nodeData (g, tf, 'moleculeType') <- 'TF'
}
cw <- CytoscapeWindow (target_cond, graph=g, overwrite=TRUE)
displayGraph (cw)
setVisualStyle(cw, 'Curved')
attribute.values <- c ('Kinase',  'TF')
node.colors      <- c ('#FEC445', '#FF3300')
setNodeColorRule (cw,
node.attribute.name='moleculeType',
control.points = attribute.values,
colors = node.colors,
mode = "lookup")
setNodeColorRule (cw,
node.attribute.name='moleculeType',
control.points = attribute.values,
colors = node.colors,
mode = "lookup")
displayGraph (cw)
setNodeColorRule (cw,
node.attribute.name='moleculeType',
control.points = attribute.values,
colors = node.colors,
mode = "lookup")
displayGraph (cw)
?setVisualStyle
?layoutNetwork
getLayoutNames()
getLayoutNames(cw)
cw <- CytoscapeWindow (target_cond, graph=g, overwrite=TRUE)
displayGraph (cw)
layoutNetwork(cw, 'circular')
setVisualStyle(cw, 'Curved')
