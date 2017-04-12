library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(xlsx)
library(dplyr)
library(growthcurver)

process_plate = function(path, meta){
    dat = read.xlsx(path,1,startRow=49,endRow=61)
    colnames(dat)[1:3] = c("Cycle","Time","Temp")
    blank = apply(dat[meta$Well[meta$Strain == 'Blank']],1,mean)
    dat$blank = blank
    dat = dat[,-c(1,3)] # Dont use them
    dat = dat[,-which(colnames(dat) %in% meta$Well[meta$Strain == 'Blank'])]
    dat$Time = dat$Time / (60*60)
    return(dat)
}

get_temporal = function(od){
    time = od[,'Time']
    blank = od$blank
    od = od[,-1]
    od = od[,-dim(od)[2]]
    od = od - blank
    od[od<=0] = min(od[od>0])
    #od = od[,-Time]
    rate = 4*apply(log2(apply(od,2,smooth)),2,diff)
    rate = rbind(rate, rep(NA, dim(rate)[2]))
    rate = as.data.frame(rate)
    # Annotate
    od$Time = time
    rate$Time = time
    od_tidy = od %>% gather(Well, OD, -Time)
    rate_tidy = rate %>% gather(Well, Growth_Rate, -Time)
    all = left_join(od_tidy,rate_tidy)
    #all = left_join(all,meta)
    return(all)
}

join_plates = function(left,right,meta){
    lm = left_join(meta,left)
    rm = left_join(meta,right)
    rm$Repeat = rm$Repeat + max(lm$Repeat)
    all  = rbind(lm, rm)
    return(all)
}
