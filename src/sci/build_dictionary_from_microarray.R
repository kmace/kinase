read_pcl = function(file) {
  return(read.table(file,
         sep='\t',
         comment.char='',
         quote="",
         header=T)[-1,-c(2,3)])
}
files = list.files('../../input/reference/spell', recursive=T, pattern='pcl',full.names =T)

cols = lapply(files, function(x) {browser(); dim(read_pcl(x))[2]})
