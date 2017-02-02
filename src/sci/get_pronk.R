library(Biobase)
library(GEOquery)
gset <- getGEO("GSE11452", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL90", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
meta = pData(gset)
px = exprs(gset)
p2t = read.table('../../input/reference/GPL_to_ORF.tsv', fill=T, header=T, stringsAsFactors=F)
rownames(px) = p2t[,2]
