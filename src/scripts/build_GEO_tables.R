rm(list=ls())
target_path = '../../intermediate/GEO/'
dir.create(file.path(target_path), showWarnings = FALSE)
source('../utils/load_libraries.R')
source('../utils/load_functions.R')
library(magrittr)
library(tidyverse)
t2g = readr::read_csv('../../intermediate/t2g.csv')#load_transcripts_to_genes()
meta = get_sample_meta()
condition_info = read_csv('../../input/meta/Condition_Meta.csv') %>% select(High_Dose, Stress) %>% unique()

meta$Strain = relevel(factor(meta$Strain), 'WT')
meta$Drug = relevel(factor(meta$Drug), 'None')
meta$Stress = relevel(factor(meta$Stress), 'None')
meta$Stress = fct_recode(meta$Stress, Glucose_Dropout = "Glucose Depletion")
meta$Media = relevel(factor(meta$Media), 'YPD')
meta = meta %>% mutate(Condition = if_else(Stress=='None',
                                           as.character(Media),
                                           as.character(Stress)
))
meta$Condition = relevel(factor(meta$Condition), 'YPD')
meta

meta = as_tibble(meta)
meta %<>% filter(Experimenter == 'Kieran')
# Build samples

meta %>%
  left_join(condition_info) %>%
  arrange(Strain, Media, Stress) %>%
  transmute(`Sample name` = Sample_Name,
         title = str_c(Sample_Name, Condition, Strain_Code, Media, Drug, sep = ':'),
         `source name` = 'Bulk population',
         organism = 'yeast - Saccharomyces cerevisiae',
         `characteristics: Strain` = str_c('w303a_', Strain, if_else(Strain == 'WT', '', '-as')),
         `characteristics: Media` = Media,
         `characteristics: Applied Stress` = Stress,
         `characteristics: Applied Stress level` = High_Dose,
         `characteristics: Inhibitors Added` = Drug,
         molecule = 'mRNA',
         `processed data file` = str_c('counts/', Sample_Name, '.tab'),
         `raw data file` = str_c('fastq/', Sample_Name, '.fastq.gz')) -> sample_table

write_csv(sample_table, str_c(target_path, 'samples.csv'))

counts_checksum = read_table('counts_checksum.tsv', col_names = FALSE)
counts_checksum$ft = '.tab'
colnames(counts_checksum) = c('file checksum', 'file name', 'file type')
counts_checksum %<>% select(`file name`,	`file type`, `file checksum`)
write_csv(counts_checksum, str_c(target_path, 'counts_checksums.csv'))

fastq_checksum = read_table('fastq_checksum.tsv', col_names = FALSE)
fastq_checksum$ft = '.fastq.gz'
colnames(fastq_checksum) = c('file checksum', 'file name', 'file type')
fastq_checksum %<>% select(`file name`,	`file type`, `file checksum`)
write_csv(fastq_checksum, str_c(target_path, 'fastq_checksums.csv'))
