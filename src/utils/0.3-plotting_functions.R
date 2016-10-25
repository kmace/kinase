# ---- _0.2 ----
library(reshape2)
library(ggplot2)

plotHclustColors <- function(matrix,labels,hang=.1,...) {
  colnames(matrix) <- labels
  d <- dist(t(matrix))
  hc <- hclust(d)
  labelColors <- brewer.pal(nlevels(labels), 'Set3')
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[which(levels(lab) == a$label)]
      attr(n, 'nodePar') <- c(a$nodePar, lab.col=labCol, pch=NA)
    }
    n
  }
  clusDendro <- dendrapply(as.dendrogram(hc,hang=hang), colLab)
  plot(clusDendro,...)
}

get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}

reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
        dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd, method = 'single')
    cormat <-cormat[hc$order, hc$order]
}

correlation_plot = function(data, upper=TRUE) {
    cormat <- cor(data)
    # Reorder the correlation matrix
    cormat <- reorder_cormat(cormat)
    upper_tri <- get_upper_tri(cormat)
    if(upper) {
        # Melt the correlation matrix
        melted_cormat <- melt(upper_tri, na.rm = TRUE)
        } else {
            # Use this if you want the full matrix
            melted_cormat <- melt(cormat, na.rm = TRUE)
        }

    # Create a ggheatmap
    g = ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
            geom_tile(color = "white") +
            scale_fill_gradient2(low = "blue",
                                 high = "red",
                                 mid = "white",
                                 midpoint = 0,
                                 limit = c(-1,1),
                                 space = "Lab",
                                 name="Pearson\nCorrelation") +
            theme_minimal() + # minimal theme
            theme(axis.text.x = element_text(angle = 45,
                                             vjust = 1,
                                             size = 7,
                                             hjust = 1)) +
            coord_fixed()
    if(dim(cormat)[1] < 100) {
        g = g +
            geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 1.5) +
            theme(axis.title.x = element_blank(),
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
    }
    return(g)
}

pca_plot = function(data, data_point_labels, label_vectors=FALSE) {
    data.pca <- prcomp(data,
                     center = TRUE,
                     scale. = TRUE)
    ggbiplot(data.pca,
             obs.scale = 1,
             var.scale = 1,
             var.axes = label_vectors,
             groups = data_point_labels,
             ellipse = TRUE,
             circle = TRUE) +
        scale_color_discrete(name = '') +
        theme(legend.direction = 'horizontal', legend.position = 'top')
}

plot_tf_expressions = function(expression, reg, file_name = 'Experssion_by_TF.pdf'){
    tfs = unique(reg$TF)
    pdf(file_name);
    for (i in 1:length(tfs)) {
        dat = na.omit(expression[reg[reg$TF == tfs[i],]$expression_index,]);
        if(length(dim(dat))>1 && dim(dat)[1]>1){
            dat[dat > 2] = 2
            dat[dat < -2] = -2
            heatmap.2(dat,
                      Colv=F,
                      main=paste(tfs[i],': (',dim(dat)[1],' genes)',sep=''),
                      col=viridis(20),
                      cexCol=.5,
                      cexRow=0.01,
                      trace='none')
        }
    }
    dev.off()
}
