source('sleuth.R')

mydata = sleuth_to_matrix(so,"obs_norm","est_counts")$data





pinc = sleuth_to_matrix(so, 'obs_norm', 'est_counts')$data
colnames(pinc)
pinc_log2 = apply(pinc,2,function(x) log2((x+1)/(1 + pinc[,"WT_NO_NO"])))
head(pinc_log2)
ls()
gasch = read.table('../Gasch/gasch-analysis/complete_dataset2.txt',header=T)
gasch = read.table('../Gasch/gasch-analysis/complete_dataset2.txt',header=T,quote='',sep='\t')
head(gasch)
ls
ls()
head(pinc)
head(gasch)
gasch %>% separate(NAME,c('Strain','Drug','Condition'),sep='_',remove=FALSE)
gasch[1,1]
gasch[1,2]
colnames(gasch[,1:3])
head(pinc_log2)
tail(pinc_log2)
pinc_log2$transcript = rownames(pinc_log2)
head(pinc_log2)
pinc_log2 = apply(pinc,2,function(x) log2((x+1)/(1 + pinc[,"WT_NO_NO"])))
pinc_log2[,'transcript'] = rownames(pinc_log2)
cbind(pinc_log2, rownames(pinc_log2))
pinc_log2 = cbind(pinc_log2, rownames(pinc_log2))
colnames(pinc_log2)
colnames(pinc_log2)[38]
colnames(pinc_log2)[38]='transcript'
colnames(pinc_log2)=="WT_NO_NO"
which(colnames(pinc_log2)=="WT_NO_NO")
pinc_log2 = pinc_log2[,-which(colnames(pinc_log2)=="WT_NO_NO")]
head(gasch)
merge(pinc_log2,gasch, by = c("transcript","UID")
)
colnames(gasch)
colnames(pinc_log2)
merge(pinc_log2,gasch, by = c("transcript","UID"))
merge(pinc_log2,gasch, by = c(row.names,"UID"))
merge(pinc_log2,gasch, by = c("transcript","UID"))
head(gasch[,"UID"])
head(pinc_log2[,"transcript"])
gasch[,"UID"] = as.character(gasch[,"UID"])
merge(pinc_log2,gasch, by = c("transcript","UID"))
merge(gasch,pinc_log2, by = c("UID","transcript"))
colnames(gasch)[1]
colnames(gasch)[1] = "transcript"


all = merge(gasch,pinc_log2, by = "transcript")

mydata = all[,4:dim(all)[2]]
cormat <- round(cor(mydata),2)

mydata = mydata %>% na.omit()

dim(mydata)
#mydata = sleuth_to_matrix(so,"obs_norm","tpm")
cormat <- round(cor(mydata),2)

melted_cormat <- melt(cormat)
#ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
#  geom_tile()
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix

melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
# ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
#  geom_tile(color = "white") +
#  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#    midpoint = 0, limit = c(-1,1), space = "Lab",
#    name="Pearson\nCorrelation") +
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 45,
#                                   vjust = 1,
#                                   size = 12,
#                                   hjust = 1))+
#  coord_fixed()
reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)


pdf('gasch_and_pincus.pdf')
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue",
                      high = "red",
                      mid = "white",
                      midpoint = 0,
                      limit = c(-1,1),
                      space = "Lab",
                      name="Pearson\nCorrelation") +
  theme_minimal() + # minimal theme
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 1,
                                   size = 2,
                                   hjust = 1),
       axis.text.y = element_text(angle = 0,
                                        vjust = 1,
                                        size = 2,
                                        hjust = 1)
                                   ) +
  coord_fixed()
ggheatmap +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.6, 0.7),
  legend.direction = "horizontal") +
guides(fill = guide_colorbar(barwidth = 7,
                             barheight = 1,
                             title.position = "top",
                             title.hjust = 0.5))
dev.off()
