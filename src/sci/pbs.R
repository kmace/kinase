pdf('pbs2.pdf')
for(s in levels(meta$Condition)) {
  plot = ggplot(meta_pc, aes(x=PC1, y=PC2)) +
    geom_point(data = subset(meta_pc, Condition == s),
               aes(fill = Condition,
                   size = if_else(Strain %in% c('TPK123','PBS2'),4,2),
                   color = if_else(Strain %in% c('TPK123','PBS2'),2,1)),
               shape = 21,
               stroke=1) +
    geom_text(aes(label=Strain),
              size=1) +
    scale_fill_discrete(drop = FALSE) +
    labs(title = s) +
    guides(size=FALSE, color=FALSE, fill=FALSE) +
    scale_size(range = c(3, 7))
  print(plot)
}
dev.off()