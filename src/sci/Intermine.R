library(InterMineR)
library(tidyverse)
load('../../intermediate/images/paper_data.RData')
im.yeast = initInterMine(listMines()["YeastMine"], token='q1w651o2B4GeG1VaT0cd')

hsa_model = getModel(im.yeast)

yeast.widgets = as.data.frame(getWidgets(im.yeast))




std_residuals %>% 
filter(abs(residual)>2.5) %>% 
mutate(sign = sign(residual)) %>% 
group_by(sample_id, sign) %>% 
left_join(t2g) %>% nest() -> grouped

grouped %>% sample_n(5) %>%
mutate(goi_gene = data %>% map('Gene'), 
       go_res = map_df(goi_gene, function(x) try(doEnrichment(im = im.yeast, 
                                                           ids = x,
                                                           widget = "go_enrichment_for_gene")$data))) -> go_residuals




goi = names(which(std_resid_matrix[,'PBS2_Salt']>2.5))
goi_gene = t2g %>% filter(name %in% goi) %>% pull(Gene)

GO_enrichResult = doEnrichment(
  im = im.yeast,
  ids = goi_gene,
  widget = "go_enrichment_for_gene"
  )

GO_enrichResult = doEnrichment(
  im = im.yeast,
  genelist = "IRE1 Tunicamycin overestimated proper",
  widget = "go_enrichment_for_gene"
  )