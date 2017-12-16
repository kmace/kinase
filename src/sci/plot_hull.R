plot_hull = function(kinase_x,kinase_y){
if(kinase_x==kinase_y){
  return(NA)
}
top_residuals %>%
    select(name, residual, Strain_Code, Condition) %>%
    group_by(name,Condition) %>%
    spread(key = Strain_Code, value = residual, fill = 0) %>%
    ungroup()-> dat2

colnames(dat2)[which(colnames(dat2)==kinase_x)] = 'x'
colnames(dat2)[which(colnames(dat2)==kinase_y)] = 'y'

s2 = dat2 %>% split(dat2$Condition)
ch2 = s2 %>%  lapply(., function(el) chull(el$x,el$y))
ch2 = lapply(names(ch2), function(x) s2[[x]][ch2[[x]],]) %>%
do.call(rbind, .)
ggplot(ch2, aes(x=x,y=y, color=Condition, group=Condition, fill=Condition)) + geom_polygon(alpha=.2) + xlab(kinase_x) + ylab(kinase_y)

}
