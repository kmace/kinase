V = na.omit(unique(c(dat$module,dat$Strain_Code)))

edLM = dat %>%
  filter(de>1.5)  %>%
  split(.$Strain_Code) %>%
  map(~ .x %>% select(-Strain_Code) %>%
        rename(edges = module, weights = de) %>%
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
