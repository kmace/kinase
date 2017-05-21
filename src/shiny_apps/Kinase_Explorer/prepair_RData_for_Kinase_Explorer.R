library(dplyr)
library(tidyr)

# Load YesTFasco Data
#load('../../input/reference/20120129_allMotifData1.02.rdat')
tf_meta = data.frame(cname = colnames(dataMat))  %>% 
  separate(cname, c('TF', 'experiment_code'), sep='_', remove=F)

tf_meta$expert = FALSE
tf_meta$expert[expert] = TRUE
tf_meta$name = t2g[match(tf_meta$TF, t2g$target_id),]$name
