read_pcl = function(file) {
    raw = read.table(file,
         sep='\t',
         comment.char='',
         quote="",
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

library(parallel)

files = list.files('../../input/reference/spell', recursive=T, pattern='pcl',full.names =T)

# Calculate the number of cores
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl,list("read_pcl"))
junk <- clusterEvalQ(cl, library(dplyr))

dims = parLapply(cl, files,
          function(x){
            dim(read_pcl(x))
            })



dat = do.call(rbind, dims)
dat = as.data.frame(dat)
dat$files = files
colnames(dat) = c('nrows','ncols','path')
