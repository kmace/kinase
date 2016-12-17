source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')

data = assay(vsd)
data = rename_gene_names(data,t2g)
col_metadata = colData(vsd)[,1:11]
row_metadata = t2g[match(rownames(data),t2g$name),]

#Filter Columns
with_drug = col_metadata$Drug == 'Cocktail'
data = data[,with_drug]
col_metadata = col_metadata[with_drug,]
colnames(data) = paste(col_metadata$Stress, col_metadata$Strain_Code, col_metadata$Media, sep='_')

#Sort Columns
stress_order = order(col_metadata$Stress)
col_metadata = col_metadata[stress_order,]
data = data[,stress_order]

gene_median = apply(data,1,median)
data = apply(data, 2, function(x) x - gene_median)
row_metadata$gene_median = gene_median

not_low_genes = which(gene_median > 2.1)
data = data[not_low_genes,]
row_metadata = row_metadata[not_low_genes,]



dir.create('../../input/gct')

writeGCT(data,
         row_metadata,
         col_metadata,
         '../../input/gct/with_drug.gct')
