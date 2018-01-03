library(RCyjs)
library(graph)
genes %>% select(data, data_augment, name) %>% unnest() %>% arrange(desc(abs(.resid))) %>% select(name, Condition, Strain_Code, Strain, Expression, .fitted, .resid, .std.resid) %>% group_by(Condition, name) %>% mutate(de = Expression - mean(Expression[Strain == 'WT'])) -> links
links %>% dplyr::filter(abs(de) > 4) %>% filter(name != 'HIS3') -> links_strong
links_strong %<>% ungroup() %>% select(Strain, name, Condition, de, .resid, .std.resid)

g = new("graphNEL", edgemode = "directed")
rcy <- RCyjs(portRange=9047:9067, quiet=TRUE, graph=g)

for(cond in unique(links_strong$Condition)[1:3]){
link_sub = links_strong[links_strong$Condition == cond, ]

g = new("graphNEL", edgemode = "directed")
nodeDataDefaults(g, attr = "type") <- "gene"

edgeDataDefaults(g, attr = "Condition") <- "none"
edgeDataDefaults(g, attr = "de") <- 0
edgeDataDefaults(g, attr = "resid") <- 0

for(kinase in unique(link_sub$Strain)){
  g = graph::addNode(kinase, g)
  nodeData(g, kinase, "type") = "kinase"
}

for(gene in unique(link_sub$name)){
  g = graph::addNode(gene, g)
}

g = graph::addEdge(as.character(link_sub$Strain), link_sub$name, g, link_sub$de)
addGraph(rcy, g)
}
redraw(rcy)
layout(rcy, "cose")
