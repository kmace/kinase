stat_smooth_func <- function(mapping = NULL, data = NULL,
                             geom = "smooth", position = "identity",
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc <- ggproto("StatSmooth", Stat,

                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))

                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }

                            params
                          },

                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }

                            if (is.null(data$weight)) data$weight <- 1

                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }

                            if (is.character(method)) method <- match.fun(method)

                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))

                            m = model
                            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                                             list(a = format(coef(m)[1], digits = 3),
                                                  b = format(coef(m)[2], digits = 3),
                                                  r2 = format(summary(m)$r.squared, digits = 3)))
                            func_string = as.character(as.expression(eq))

                            if(is.null(xpos)) xpos = min(data$x)*0.9
                            if(is.null(ypos)) ypos = max(data$y)*0.9
                            data.frame(x=xpos, y=ypos, label=func_string)

                          },

                          required_aes = c("x", "y")
)

per_condition_subset_results = read_csv('../../intermediate/per_condition_subset_results.csv')
wt_results = read_csv('../../intermediate/wt_results.csv')


salt_genes = wt_results %>% dplyr::filter(condition=='Salt') %>% pull(name)
per_condition_subset_results %>%
  dplyr::filter(condition == 'Salt', name != 'HIS3') %>%
  dplyr::filter(name %in% salt_genes) %>%group_by(name) %>% mutate(n=n()) %>% arrange(desc(n)) %>% dplyr::filter(n<10) %>%
  select(Kinase, name) %>%
  split(.$Kinase) %>%
  lapply(function(x) pull(x,name)) -> sets

library(SuperExactTest)
res = supertest(sets, n = length(salt_genes), degree = c(2:4))
tab = summary(res)$Table %>% as_tibble() %>% arrange(P.value)  %>% mutate(log2FE = log2(FE))
tab %>% arrange(P.value)  %>% mutate(log2FE = log2(FE)) %>% `[`(1,7) %>% pbcopy()

good = res$P.value < 0.01 & res$overlap.sizes > 10
res$P.value = res$P.value[good]
res$overlap.sizes = res$overlap.sizes[good]
plot(res,degree=2, sort.by='p-value')


statistic = conditional_results %>% dplyr::filter(condition == 'Salt', Kinase == 'PBS2')
module.set = modules %>% split(.$module) %>% lapply(function(x) pull(x, name))
idx = ids2indices(module.set, statistic$name, remove.empty=TRUE)
cameraPR(statistic$stat, idx) %>% as_tibble(rownames = 'module') -> out
out
#spread(key=Kinase, value = pvalue) %>%
#as.data.frame() %>% remove_rownames() %>%
#column_to_rownames('name') %>% is.na() %>% `!` -> mat

# by Condition
cond = 'Tunicamycin'

per_condition_subset_results %>%
  left_join(wt_results %>%
              dplyr::select(condition, log2FoldChange, name, padj) %>%
              dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
  dplyr::filter(condition == cond & wtp<0.1) %>%
  ggplot(aes(x=wt_change, y=log2FoldChange)) +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point() + geom_smooth(method = 'lm') + facet_wrap(~Kinase)

for(cond in levels(meta$Condition)[-1]){
  pdf(paste0(cond,'.pdf'))
  p = per_condition_subset_results %>%
    left_join(wt_results %>%
                dplyr::select(condition, log2FoldChange, name, padj) %>%
                dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
    dplyr::filter(condition == cond & wtp<0.1) %>%
    mutate(side = sign(wt_change)) %>%
    ggplot(aes(x=wt_change, y=log2FoldChange, group=side)) +
    stat_poly_eq(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "\n")),
                 formula = formula, parse = TRUE, coef.digits = 2, label.y.npc = "bottom") +
    geom_point() + geom_smooth(method='lm') + facet_wrap(~Kinase)
  print(p)
  dev.off()
}


# by Strain
kinase = 'CDC15'

per_condition_subset_results %>%
  left_join(wt_results %>%
              dplyr::select(condition, log2FoldChange, name, padj) %>%
              dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
  dplyr::filter(Kinase == kinase & wtp<0.1 & abs(wt_change)>.5)%>%
  ggplot(aes(x=wt_change, y=log2FoldChange)) +
  stat_smooth_func(geom="text",hjust=0,parse=TRUE) +
  geom_point() + geom_smooth() + facet_wrap(~condition)


# by Strain
kinase = 'HOG1'
for(kinase in levels(meta$Strain)[-1]){
pdf(paste0(kinase,'.pdf'))
p = per_condition_subset_results %>%
  left_join(wt_results %>%
              dplyr::select(condition, log2FoldChange, name, padj) %>%
              dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
  dplyr::filter(Kinase == kinase & wtp<0.1)%>%
  mutate(side = sign(wt_change)) %>%
  ggplot(aes(x=wt_change, y=log2FoldChange, group=side)) +
  stat_poly_eq(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "\n")),
                formula = formula, parse = TRUE, coef.digits = 2, label.y.npc = "bottom") +
  geom_point() + geom_smooth(method='lm') + facet_wrap(~condition)
print(p)
  dev.off()
}


genes_de_once = wt_results %>% filter(padj<0.01 & abs(log2FoldChange)>.5) %>% pull(name) %>% unique()
sample_meta = colData(dds_wt)
rld = rlog(dds_wt)
rlog = assay(rld)
rlog = rlog[genes_de_once, ]
source('../paper/figures/colors.R')


library(ComplexHeatmap)

sample_ha =  HeatmapAnnotation(sample_meta %>% as.data.frame %>% select(Condition),
                               col = master_col,
                               show_annotation_name = TRUE)

Heatmap(rlog - rowMeans(rlog[,sample_meta$Condition == 'YPD']),
        top_annotation = sample_ha,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        use_raster = F,
        show_column_names = F)

x = selectArea()
rownames(rlog)[x$row_order] %>% pbcopy


test = coef(dds_wt)
test = test[genes_de_once,]
colnames(test) %<>% word(2,sep = '_')
test = test[,-1]
Heatmap(test,
        #km = 5,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        use_raster = F)

x = selectArea()
rownames(test)[x$row_order] %>% pbcopy


