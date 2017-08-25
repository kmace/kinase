source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')
x = assay(vst)
#x = x - rowMeans(x)


gasch = load_gasch()
gx = gasch$data
gmeta = gasch$meta

pronk = load_pronk()
px = pronk$data
pmeta = pronk$meta

mega = load_megaYeast()
mx = mega$data
mmeta = mega$meta

common_genes = Reduce(intersect,
                      list(rownames(gx),
                           rownames(x),
                           rownames(px),
                           rownames(mx)))

all_genes = Reduce(union,
                      list(rownames(gx),
                           rownames(x),
                           rownames(px),
                           rownames(mx)))

venn(list(gasch = rownames(gx),
          mydata = rownames(x),
          pronk = rownames(px),
          mega = rownames(mx)))

'%!in%' <- function(x,y)!('%in%'(x,y))

missing = venn(list(gasch = all_genes[all_genes %!in% rownames(gx)],
          mydata = all_genes[all_genes %!in% rownames(x)],
          pronk = all_genes[all_genes %!in% rownames(px)],
          mega = all_genes[all_genes %!in% rownames(mx)]))
intersections = attr(missing,'intersections')

x = x[rownames(x) %!in% intersections$`gasch:pronk:mega`, ]

datasetlist = list(
  scale(t(x)),
  scale(t(px)),
  scale(t(gx)),
  scale(t(mx))
)

library(INSPIRE)

res = INSPIRE(datasetlist = datasetlist, mcnt = 100, lambda = 0.01, printoutput = 1)


clusters = res$Z[[1]]
clusters = tibble(Gene = names(clusters), Cluster = clusters)
clusters = left_join(clusters, t2g, by=c('Gene' = 'target_id'))
clusters = clusters %>% arrange(Cluster)

dir.create('../../input/clusters/INSPIRE_clusters')
clusters %>%
  group_by(Cluster) %>%
  do(write_csv(., path = paste0('../../input/clusters/INSPIRE_clusters/cluster_',
                                as.character(unique(.$Cluster)),
                                '.csv')))

write.csv(clusters,file='../../input/clusters/INSPIRE_clusters/all_clusters.csv')

