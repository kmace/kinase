library(dplyr)
library(tidyr)
select = dplyr::select
meta = read.table('../meta/Pincus/hiseq_info.txt', header=T, stringsAsFactors=FALSE, sep = '\t')
meta = meta %>%
    separate(Barcode,c('Tnum','bar'),sep='-') %>%
    mutate(filename = paste(bar,'-s_',as.character(Position),'_1_sequence.txt.tar.gz',sep = ''),
           path = file.path("../data/Pincus_data/results",filename)) %>%
    select(sample = SampleName,
           path)  %>%
    separate(sample,c('Strain','Drug','Condition'),sep='_',remove=FALSE) %>%
    mutate(Strain = relevel(factor(Strain),'WT'),
           Drug = relevel(factor(Drug),'NO'),
           Condition = relevel(factor(Condition),'NO')) %>%
    filter(Condition=='NO', Drug=='NO')

n_experiments = dim(meta)[1]


qval = list()


for (i in 1:200) {
random_perms = matrix(sample(c(TRUE,FALSE), n_experiments*1, replace = TRUE),n_experiments,1)
#random_perms = matrix(sample(c(TRUE,FALSE), n_experiments*200, replace = TRUE),n_experiments,200)
if(sum(random_perms)>1 && sum(random_perms)<n_experiments) {
model = model.matrix(~., data=data.frame(random_perms))
so <- sleuth_prep(meta, model)
#for (i in 1:200) {
so <- sleuth_fit(so)
so = sleuth_wt(so, which_beta = 'random_permsTRUE')
qval = append(qval,(sleuth_results(so,'random_permsTRUE') %>% arrange(qval) %>% select(qval))[[1]][1])
}}
