source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')

# only look at WT Cells
wt = rlog_meta$Strain == 'WT'
rlog = rlog[ ,wt]
rlog_meta = rlog_meta[wt, ]

# only look at YPD (we only have no drug control for YPD)
ypd = rlog_meta$Media == 'YPD'
rlog = rlog[ ,ypd]
rlog_meta = rlog_meta[ypd, ]

get_sample_mean = function(x) apply(x,1,mean)
# Establish Control sets, these will be denominators throughout.
control_noDrug = which(
    rlog_meta$Stress == 'None' &
    rlog_meta$Drug == 'None')

control_Drug = which(
    rlog_meta$Stress == 'None' &
    rlog_meta$Drug == 'Cocktail')

control_noDrug = get_sample_mean(rlog[,control_noDrug])
control_Drug = get_sample_mean(rlog[,control_Drug])

availability_table = data.frame(rlog_meta) %>% dplyr::select(Stress, Drug) %>% table
have_data = apply(availability_table,1,min) != 0
conditions = rownames(availability_table)[have_data]
# now for heatshock
plot_comp = function(condition = 'Heatshock') {
    condition_noDrug = which(
        rlog_meta$Stress == condition &
        rlog_meta$Drug == 'None')

    condition_Drug = which(
        rlog_meta$Stress == condition &
        rlog_meta$Drug == 'Cocktail')

    condition_noDrug = get_sample_mean(rlog[,condition_noDrug])
    condition_Drug = get_sample_mean(rlog[,condition_Drug])

    drug = condition_Drug - control_Drug
    noDrug = condition_noDrug - control_noDrug

    c = cor(drug, noDrug)

    plot(noDrug, # X
         drug, # Y
         xlab='FC without Drug',
         ylab='FC with Drug',
         main = paste('Comparing ', condition, ' vs. YPD: with and without Drug. cor = ', c)
         )
    abline(0,1)
}

dir.create('../../output/compare_drug_to_no_drug.R')
pdf('../../output/compare_drug_to_no_drug.R/drug_normalization.pdf')
lapply(conditions, plot_comp)
dev.off()



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
