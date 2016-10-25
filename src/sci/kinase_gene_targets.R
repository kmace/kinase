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
           Condition = relevel(factor(Condition),'NO')) #%>%
    #filter(Condition=='HS')

library("rhdf5")
library("sleuth")
library("biomaRt")

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl")
mart
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)


t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

model = model.matrix(~ 1 + Strain, data = meta)
so <- sleuth_prep(meta, model, target_mapping = t2g)
so <- sleuth_fit(so)

so <- sleuth_wt(so, which_beta = 'StrainATG1')
so <- sleuth_wt(so, which_beta = 'StrainFUS3')
so <- sleuth_wt(so, which_beta = 'StrainHOG1')
so <- sleuth_wt(so, which_beta = 'StrainIRE1')
so <- sleuth_wt(so, which_beta = 'StrainKSS1')
so <- sleuth_wt(so, which_beta = 'StrainPKC1')
so <- sleuth_wt(so, which_beta = 'StrainRIM15')
so <- sleuth_wt(so, which_beta = 'StrainSCH9')
so <- sleuth_wt(so, which_beta = 'StrainSLT2')
so <- sleuth_wt(so, which_beta = 'StrainSNF1')
so <- sleuth_wt(so, which_beta = 'StrainTPK')
