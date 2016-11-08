source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')


# TODO create overlay for genes of interest

plot_drug_comp = function(condition = 'Heatshock') {
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

    corr = cor(drug, noDrug)
    #covv = cov(drug, noDrug)
    plot(noDrug, # X
         drug, # Y
         xlab='FC without Drug',
         ylab='FC with Drug',
         main = paste('The FC between ', condition,
                      ' and YPD.\n Comparing with and without Drug.\n',
                      'cor = ', round(corr,2))#,
                      #'\ncov = ', round(covv,3))
         )
    abline(0,1)
    abline(v=0,lty=3)
    abline(h=0,lty=3)
    d = drug[t2g$target_id[t2g$name %in% hs]]
    nd = noDrug[t2g$target_id[t2g$name %in% hs]]
    points(nd,d, col='red', pch=16)

}

plot_cond_comp = function(condition = 'Heatshock') {
    condition_noDrug = which(
        rlog_meta$Stress == condition &
        rlog_meta$Drug == 'None')

    condition_Drug = which(
        rlog_meta$Stress == condition &
        rlog_meta$Drug == 'Cocktail')

    condition_noDrug = get_sample_mean(rlog[,condition_noDrug])
    condition_Drug = get_sample_mean(rlog[,condition_Drug])

    condition_val = condition_Drug - condition_noDrug
    control_val =  control_Drug - control_noDrug

    corr = cor(condition_val, control_val)
    covv = cov(condition_val, control_val)
    plot(control_val, # X
         condition_val, # Y
         xlab='Drug FC in YPD',
         ylab=paste('Drug FC in ', condition),
         main = paste('The FC between Drug and no Drug.\n',
                      'Comparing: with and without ', condition, '.\n',
                      'cor = ', round(corr,2))#,
                      #'cov = ', round(covv,3))
         )
    abline(0,1)
}

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

hs = load_gene_list_from_file('../../input/reference/HSF1_targets_in_HS_yestract.txt')

dir.create('../../output/compare_drug_to_no_drug.R')
pdf('../../output/compare_drug_to_no_drug.R/drug_vs_noDurg.pdf')
lapply(conditions, plot_drug_comp)
dev.off()
pdf('../../output/compare_drug_to_no_drug.R/Stress_vs_noStress.pdf')
lapply(conditions, plot_cond_comp)
dev.off()



# hs = apply(assay(rld_all[HS_target_id,hs]),1,mean)
#
# t = cbind(control, hs)
# colnames(t) = c('c','c','c','c','hs','hs','hs','hs')
# tt = t(apply(t,1,function(x) x - mean(x)))
# pheatmap(tt, main='None')
#
# dev.new()
# control = which(
#     colData(rld_all)$Stress == 'None' &
#     colData(rld_all)$Drug == 'Cocktail' &
#     colData(rld_all)$Media == 'YPD' &
#     colData(rld_all)$Strain == 'WT')
#
# hs = which(
#     colData(rld_all)$Stress == 'Heatshock' &
#     colData(rld_all)$Drug == 'Cocktail' &
#     colData(rld_all)$Media == 'YPD' &
#     colData(rld_all)$Strain == 'WT')
#
# control = assay(rld_all[HS_target_id,control])
# hs = assay(rld_all[HS_target_id,hs])
