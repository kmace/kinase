source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')

control_noDrug = which(
    rlog_meta$Stress == 'None' &
    rlog_meta$Drug == 'None' &
    rlog_meta$Media == 'YPD' &
    rlog_meta$Strain == 'WT')

condition_noDrug = which(
    rlog_meta$Stress == 'Heatshock' &
    rlog_meta$Drug == 'None' &
    rlog_meta$Media == 'YPD' &
    rlog_meta$Strain == 'WT')

control_Drug = which(
    rlog_meta$Stress == 'None' &
    rlog_meta$Drug == 'Cocktail' &
    rlog_meta$Media == 'YPD' &
    rlog_meta$Strain == 'WT')

condition_Drug = which(
    rlog_meta$Stress == 'Heatshock' &
    rlog_meta$Drug == 'Cocktail' &
    rlog_meta$Media == 'YPD' &
    rlog_meta$Strain == 'WT')

get_sample_mean = function(x) apply(x,1,mean)

control_noDrug = get_sample_mean(rlog[,control_noDrug])
condition_noDrug = get_sample_mean(rlog[,condition_noDrug])
control_Drug = get_sample_mean(rlog[,control_Drug])
condition_Drug = get_sample_mean(rlog[,condition_Drug])

c = cor(condition_Drug   - control_Drug,
        condition_noDrug - control_noDrug)

plot(condition_Drug   - control_Drug,
     condition_noDrug - control_noDrug)




hs = apply(assay(rld_all[HS_target_id,hs]),1,mean)

t = cbind(control, hs)
colnames(t) = c('c','c','c','c','hs','hs','hs','hs')
tt = t(apply(t,1,function(x) x - mean(x)))
pheatmap(tt, main='None')

dev.new()
control = which(
    colData(rld_all)$Stress == 'None' &
    colData(rld_all)$Drug == 'Cocktail' &
    colData(rld_all)$Media == 'YPD' &
    colData(rld_all)$Strain == 'WT')

hs = which(
    colData(rld_all)$Stress == 'Heatshock' &
    colData(rld_all)$Drug == 'Cocktail' &
    colData(rld_all)$Media == 'YPD' &
    colData(rld_all)$Strain == 'WT')

control = assay(rld_all[HS_target_id,control])
hs = assay(rld_all[HS_target_id,hs])
