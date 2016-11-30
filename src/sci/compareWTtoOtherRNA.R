rm(list=ls())
source('../utils/load_libraries.R')
source('../utils/load_functions.R')

t2g = load_transcripts_to_genes()

meta = get_sample_meta()
raw_counts = load_count_matrix(meta,'star')

lib_size = apply(raw_counts,2,sum)
p = t(apply(raw_counts,1,function(x) x/lib_size))
m = p / apply(p,1,function(x) sqrt(sum(x^2)))
NoDrug_joanna = h[,1]/c[,1]
Drug_26 = h[,"26h"]/c[,"26yc"]
Drug_25 = h[,"25h"]/c[,"25yd"]
Drug_26 = h[,"26h"]/c[,"26yd"]
Drug_27 = h[,"27h"]/c[,"27yd"]
Drug_28 = h[,"28h"]/c[,"28yd"]
No_Drug_25 = h[,"7k"]/c[,"25yc"]
No_Drug_26 = h[,"15k"]/c[,"26yc"]
No_Drug_27 = h[,"23k"]/c[,"27yc"]
No_Drug_28 = h[,"31k"]/c[,"28yc"]
comp = cbind(Drug_25,Drug_26,Drug_27,Drug_28,No_Drug_25,No_Drug_26,No_Drug_27,No_Drug_28,NoDrug_joanna)
pairs(comp)
pairs(log(comp))
correlation_plot(log2(comp))

cor(log2(comp),method="spear")

cormat = cor(log2(comp),method="spear")
upper=TRUE
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


pdf('Heatshock Baseline Correction.pdf')
g + ggtitle('Heatshock Data')
dev.off()
