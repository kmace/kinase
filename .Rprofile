len = length
pbcopy = function(x,sep="\t",col.names=T,...) {
  write.table(x
             ,file = pipe("pbcopy")
             ,sep=sep
             ,col.names = col.names
             ,row.names = F
             ,quote = F,...)
}

pbpaste = function(sep="\t",header=T,...) {
  read.table(pipe("pbpaste")
            ,sep=sep
            ,header=header,...)
}
