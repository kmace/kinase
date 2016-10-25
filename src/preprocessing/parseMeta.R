args = commandArgs(trailingOnly=TRUE)
xls_file = args[1]
meta = read.table(xls_file, header=T, stringsAsFactors=FALSE, sep = '\t')
library(dplyr)
library(tidyr)
meta = meta %>%
        separate(Barcode,c('Tnum','bar'),sep='-') %>%
        mutate(filename = paste(bar,'-s_',as.character(Position),'_1_sequence.txt.tar.gz',sep = '')) %>%
	select(SampleName, filename)

write.table(meta, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=',', file = './meta/meta.tsv')
