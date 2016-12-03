
library(DESeq2)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(viridis)
library(tidyr)
library(genefilter)
library(ggrepel)
library(plotly)
library(Glimma)

source('0.1-data_loading_functions.R')
source('0.2-data_manipulation_functions.R')
source('0.3-plotting_functions.R')

t2g = load_transcripts_to_genes()
sample_meta = load_sample_meta_data()
seq_meta = load_seq_meta_data()
meta = right_join(sample_meta, seq_meta, by = c("Sample_Name" = "sample"))
raw_counts = load_count_matrix(meta,'tophat')

colnames(meta)

model_formula = ~ 1
#model_formula = ~ Stress

deseq_object = load_deseq_object(model_formula, meta, raw_counts)

glMDPlot(dds)

deseq_counts = counts(deseq_object, normalized=TRUE)
colnames(deseq_counts) = meta$Sample_Name

fc = filter_and_foldchange(deseq_counts)
correlation_plot(fc)

deseq_object.vst = varianceStabilizingTransformation(deseq_object, blind=TRUE, fitType="mean")
deseq_object.vst %>% head
x = assay(deseq_object.vst)
rv = rowVars(x)
rv = genefilter::rowVars(x)
select = order(rv, decreasing = TRUE)[seq_len(6800)]
pca = prcomp(t(x[select, ]))

#head(cbind(meta$Sample_Name,rownames(pca$x)[match(meta$Sample_Name,rownames(pca$x))]))
meta_pc = cbind(meta,pca$x[match(meta$Sample_Name,rownames(pca$x)),])


library(repr)

p1 = ggplot(meta_pc, aes(x=PC1, y=PC2)) 
p2 = ggplot(meta_pc, aes(x=PC1, y=PC3)) 
p3 = ggplot(meta_pc, aes(x=PC2, y=PC3)) 
pdf('test.pdf')
options(repr.plot.width=11, repr.plot.height=10)
p1 + geom_point(shape = 21, aes(fill = Stress, colour=Drug), size=4, stroke=1) + geom_text(aes(label=Strain), size=1)
p2 + geom_point(shape = 21, aes(fill = Stress, colour=Drug), size=4, stroke=1) + geom_text(aes(label=Strain), size=1)
p3 + geom_point(shape = 21, aes(fill = Stress, colour=Drug), size=4, stroke=1) + geom_text(aes(label=Strain), size=1)
dev.off()

head(meta_pc)

p <- plot_ly(meta_pc, x = ~PC1, y = ~PC2) %>%
             add_markers(color = ~Stress, visible = T) %>%
             add_markers(color = ~Drug, visible = F) %>%
             add_markers(color = ~Experimenter, visible = F)

p <- p %>% layout(
  title = "Drop down menus - Styling",
  #xaxis = list(domain = c(0.1, 1)),
  #yaxis = list(title = "y"),
  updatemenus = list(
    list(
      y = 0.7,
      buttons = list(
        list(method = "restyle",
             args = list("visible", list(T, T, T)),
             label = "Stress"),
          
        list(method = "restyle",
             args = list("visible", list(F, T, F)),
             label = "Drug"),
 
        list(method = "restyle",
             args = list("visible", list(F, F, T)),
             label = "Experimenter")))
  ))


embed_notebook(p)

c = raw_counts[,which(meta$Stress == 'None' & meta$Drug == 'None')]
h = raw_counts[,meta$sample[which(meta$Stress == 'Heatshock' & meta$Drug == 'None')]]

library(infotheo)
library(entropy)

plot(apply(pca$x,2,var), type='l')
plot(apply(pca$x,2,function(x)  entropy(discretize(x, numBins=5, r=range(pca$x)), unit="log2")), type='l')
plot(apply(pca$x,2,function(x)  entropy(discretize(x, numBins=5), unit="log2")), type='l')

#hist(pca$x[,170])

information_plot = function (data, upper = TRUE) 
{
    cormat <- cor(data)
    cormat <- reorder_cormat(cormat)
    upper_tri <- get_upper_tri(cormat)
    if (upper) {
        melted_cormat <- melt(upper_tri, na.rm = TRUE)
    }
    else {
        melted_cormat <- melt(cormat, na.rm = TRUE)
    }
    g = ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + 
        geom_tile(color = "white") + scale_fill_gradient2(low = "blue", 
        high = "red", mid = "white", midpoint = 0, limit = c(-1, 
            1), space = "Lab", name = "Pearson\nCorrelation") + 
        theme_minimal() + theme(axis.text.x = element_text(angle = 45, 
        vjust = 1, size = 7, hjust = 1)) + coord_fixed()
    if (dim(cormat)[1] < 100) {
        g = g + geom_text(aes(Var2, Var1, label = round(value, 
            2)), color = "black", size = 1.5) + theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(), panel.grid.major = element_blank(), 
            panel.border = element_blank(), panel.background = element_blank(), 
            axis.ticks = element_blank(), legend.justification = c(1, 
                0), legend.position = c(0.6, 0.7), legend.direction = "horizontal") + 
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1, 
                title.position = "top", title.hjust = 0.5))
    }
    return(g)
}

pca$rotation %>% dim


pca$x %>% dim


sort(pca$rotation[,1], decreasing = T)

meta$Stress

renamed = deseq_counts
colnames(renamed) = unite(meta, j, Strain_Code, Plate_Code)$j
correlation_plot(filter_and_foldchange(renamed[,meta$Stress == '']))

colors = rainbow(length(unique(iris$Species)))

colors = rainbow(length(unique(iris$Species)))

library(tsne)
colors = rainbow(length(unique(meta$Plate_Code)))
names(colors) = unique(meta$Plate_Code)
ecb = function(x,y){ 
    plot(x,t='n'); 
    text(x,labels=meta$Plate_Code, col=colors[meta$Plate_Code]) 
} 

tsne_iris = tsne(t(x[select, ]), epoch_callback = ecb, perplexity=50)


