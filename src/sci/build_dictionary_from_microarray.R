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

get_pc = function(pcl, num=5) {

    if(dim(pcl)[2]<num+1){
        return(NULL)
    }
    return(tryCatch(prcomp(t(pcl))$rotation[,1:num],error=function(e) NULL))
}

get_stats = function(pcl){
  min = apply(pcl,1,function(x) min(x, na.rm=T))
  mean = apply(pcl,1,function(x) mean(x, na.rm=T))
  med = apply(pcl,1,function(x) median(x, na.rm=T))
  max = apply(pcl,1,function(x) max(x, na.rm=T))
  sd = apply(pcl,1,function(x) sd(x, na.rm=T))
  range = max - min
  return(cbind(min,mean,med,max,sd,range))
}

mean_sub = function(pcl){
  gene_average = apply(pcl,1,mean)
  pcl_mean_sub = apply(pcl,2,function(x) x - gene_average)
  return(pcl_mean_sub)
}

files = list.files('../../input/reference/spell', recursive=T, pattern='pcl',full.names =T)

# Calculate the number of cores
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl,list("read_pcl", "get_pc", "get_stats", "mean_sub"))
junk <- clusterEvalQ(cl, library(dplyr))

pcls  = parLapply(cl, files, read_pcl)
pcls = parLapply(cl, pcls, mean_sub)
stats = parLapply(cl, pcls, get_stats)
pcs   = parLapply(cl, pcls, get_pc)
files = files[!unlist(lapply(pcs,is.null))]
pcs = pcs[!unlist(lapply(pcs,is.null))]

for (i in seq(1:length(files))) {
  colnames(pcs[[i]]) = paste(files[i],1:5,sep='_')
}
pc_df_t = parLapply(cl, pcs, function(x) data.frame(t(x)))
dictionary = do.call(rbind.fill,pc_df_t)
dictionary = t(dictionary)
colnames(dictionary) = unlist(lapply(pcs, colnames))
assay_num = apply(dictionary,1,function(x) length(which(is.na(x))))
dictionary = dictionary[assay_num<5*400,]
dictionary[is.na(dictionary)] = 0

save.image('../../input/images/dictionary.RData')

# dat = do.call(rbind, dims)
# dat = as.data.frame(dat)
# dat$files = files
# colnames(dat) = c('nrows','ncols','path')
