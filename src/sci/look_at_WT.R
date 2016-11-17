source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')


# Get WTs
raw_counts = raw_counts[ ,meta$Strain=='WT']
dds = dds[,colData(dds)$Strain == 'WT']
meta = meta[meta$Strain=='WT',]
dds = estimateSizeFactors(dds)
dds_counts = counts(dds, normalized=TRUE)



conditions = unique(meta$Condition)

 plot_ll = function(counts, x, baseline){
  base = rowMeans(counts[,meta$Condition==baseline],na.rm=T)
  cond = rowMeans(counts[,meta$Condition==x],na.rm=T)
  bad = base==0 | cond==0
  base = base[!bad]
  cond = cond[!bad]
  up = which(log2(cond) - log2(base) > 1)
  down = which(log2(base) - log2(cond) > 1)
  symbol=16
  size=.3
  plot(
    base,
    cond,
    log="xy",
    xlab=baseline,
    ylab=x,
    pch=symbol,
    cex=size,
    col = 'gray'
  )
  abline(0,1)
  points(
    base[up],
    cond[up],
    pch=symbol,
    cex=size,
    col='red'
  )
  points(
    base[down],
    cond[down],
    pch=symbol,
    cex=size,
    col='blue'
  )
}


# Damn the things I am plotting are actually from a single plate,  that is
# Plate K, could it just be that we messed up plate K? that could be the reason!!!
dir.create('../../output/look_at_WT.R')
baseline = "None_YPD_None"
pdf('../../output/look_at_WT.R/Scatter_vs_Baseline.pdf')
lapply(conditions, function(x) plot_ll(dds_counts,x,baseline))
dev.off()
baseline = "None_YPD_Cocktail"
pdf('../../output/look_at_WT.R/Scatter_vs_Drug_Baseline.pdf')
lapply(conditions, function(x) plot_ll(dds_counts,x,baseline))
dev.off()
