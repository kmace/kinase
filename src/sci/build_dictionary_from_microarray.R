library(dplyr)
library(parallel)

read_pcl = function(file) {
    raw = read.table(file,
         sep='\t',
         comment.char='',
         quote="",
         na.strings = c("NA", "null"),
         fill=T,
         header=T)[-1,-c(2,3)]
    colnames(raw)[1] = 'YORF'
    if(sum(duplicated(raw$YORF))>0){
        gene_average = raw %>% group_by(YORF) %>% summarize_all(mean)
    } else {
        gene_average = raw
    }
    gene_average = as.data.frame(gene_average)
    rownames(gene_average) = gene_average$YORF
    gene_average = gene_average[,-1, drop=F]
    return(gene_average)
}

get_pc = function(pcl, num=3) {
    if(dim(pcl)[2]<num+1){
        return(NULL)
    }
    return(prcomp(t(pcl))$rotation[,1:num])
}

get_stats = function(pcl){
  min = apply(pcl,2,function(x) min(x, na.rm=T))
  mean = apply(pcl,2,function(x) mean(x, na.rm=T))
  med = apply(pcl,2,function(x) median(x, na.rm=T))
  max = apply(pcl,2,function(x) max(x, na.rm=T))
  sd = apply(pcl,2,function(x) sd(x, na.rm=T))
  return(cbind(min,mean,med,max,sd))
}

files = list.files('../../input/reference/spell', recursive=T, pattern='pcl',full.names =T)

# Calculate the number of cores
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl,list("read_pcl", "get_pc", "get_stats"))
junk <- clusterEvalQ(cl, library(dplyr))

pcls  = parLapply(cl, files, read_pcl)
dims  = parLapply(cl, pcls, dim)
stats = parLapply(cl, pcls, get_stats)


pcs   = parLapply(cl, pcls, get_pc)


dat = do.call(rbind, dims)
dat = as.data.frame(dat)
dat$files = files
colnames(dat) = c('nrows','ncols','path')
