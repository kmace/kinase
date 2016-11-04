control_noDrug = which(
    colData(rld_all)$Stress == 'None' &
    colData(rld_all)$Drug == 'None' &
    colData(rld_all)$Media == 'YPD' &
    colData(rld_all)$Strain == 'WT')

condition_noDrug = which(
    colData(rld_all)$Stress == 'Heatshock' &
    colData(rld_all)$Drug == 'None' &
    colData(rld_all)$Media == 'YPD' &
    colData(rld_all)$Strain == 'WT')

control_Drug = which(
    colData(rld_all)$Stress == 'None' &
    colData(rld_all)$Drug == 'Cocktail' &
    colData(rld_all)$Media == 'YPD' &
    colData(rld_all)$Strain == 'WT')

condition_Drug = which(
    colData(rld_all)$Stress == 'Heatshock' &
    colData(rld_all)$Drug == 'Cocktail' &
    colData(rld_all)$Media == 'YPD' &
    colData(rld_all)$Strain == 'WT')

get_sample_mean = function(x) apply(x,1,mean)

control_noDrug = get_sample_mean(assay(rld_all[,control_noDrug]))
condition_noDrug = get_sample_mean(assay(rld_all[,condition_noDrug]))
control_Drug = get_sample_mean(assay(rld_all[,control_Drug]))
condition_Drug = get_sample_mean(assay(rld_all[,condition_Drug]))

plot(condition_Drug   - control_Drug,
     condition_noDrug - control_noDrug)

cor(condition_Drug   - control_Drug,
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
