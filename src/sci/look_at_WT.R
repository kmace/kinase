source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')
raw_counts = raw_counts[ ,meta$Strain=='WT']
meta = meta[meta$Strain=='WT',]
conditions = unique(meta$Condition)
#baseline = "None_YPD_None"
baseline = "None_YPD_Cocktail"

 plot_ll = function(x){
  base = rowMeans(raw_counts[,meta$Condition==baseline],na.rm=T)
  cond = rowMeans(raw_counts[,meta$Condition==x],na.rm=T)
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

pdf('test.pdf')
lapply(conditions, plot_ll)
dev.off()
