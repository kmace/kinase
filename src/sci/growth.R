library(tidyverse)
library(RColorBrewer)

meta = read.csv('../../input/meta/Growth_Rate_Experimental_Designs/Whitehead_ALL_in_ALL_Conditions_from_Saturation/384_Well_meta.csv',header=T)
od = read.table('../../input/growth/Whitehead_ALL_in_ALL_Conditions_from_Saturation/from_saturation.txt', header=T)
od = od[,-c(1,2)]
rate = 4*apply(log2(apply(od,2,smooth)),2,diff)
meta = meta %>% mutate(well_name = paste(Row, Column,sep=''))

dir.create('../../output/growth.R')

pdf('../../output/growth.R/growth_by_strain.pdf')
for(s in sort(as.character(unique(meta$Strain)))){
    strain = s;
    m = filter(meta, Strain==strain);
    matplot((0:(dim(rate)[1]-1) * 15)/60,
            apply(rate[,m$well_name],2,smooth),
            type='l',
            col=brewer.pal(8, "Set2")[factor(m$Condition)],
            lty='solid',
            ylab='Doublings per hour',
            ylim=c(-0.1,0.6),
            lwd=1.5,
            xlab='Time (h)',
            main = strain);
    legend('topright',legend=m$Condition,
           col=brewer.pal(8, "Set2")[factor(m$Condition)],
           lty='solid',
           lwd=6)
}
dev.off()


pdf('../../output/growth.R/growths_by_condition.pdf')
for(s in as.character(unique(meta$Condition))){
    condition = s;
    m = filter(meta, Condition==condition);
    matplot((0:93 * 15)/60,
            apply(rate[,m$well_name],2,smooth),
            type='l',
            col=brewer.pal(11, "Set3")[as.numeric(factor(m$Strain)) %% 11 + 1],
            lty=c('solid','dashed', 'dotted')[as.numeric(factor(m$Strain)) %% 3 + 1],
            lwd=1.5,
            ylab='Doublings per hour',
            ylim=c(-0.1,0.6),
            xlab='Time (h)',
            main = condition);
    legend('topleft',legend=m$Strain,
           col=brewer.pal(11, "Set3")[as.numeric(factor(m$Strain)) %% 11 + 1],
           lty=c('solid','dashed', 'dotted')[as.numeric(factor(m$Strain)) %% 3 + 1],
           lwd=1.5,
           cex=0.5)
}
dev.off()

ave_top_speed = apply(rate,2,function(x) mean(sort(x,decreasing=TRUE)[1:3]))
dat = cbind(meta, ave_top_speed[match(meta$well_name, names(ave_top_speed))])
colnames(dat)[dim(dat)[2]] = 'top_speed'

# Here I step in to make a supp figure for the paper:

dat = dat %>% group_by(Strain) %>% mutate(top_speed_n = top_speed/top_speed[Condition == 'Drug']) %>% ungroup()
ggplot(dat, aes(x = Strain, y = top_speed_n)) + geom_bar(stat="identity") + facet_wrap(~Condition, scales = 'free_y') + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4))
ggsave('~/Desktop/growth_like_slope.pdf')
write_csv(dat, '~/Desktop/growth_rates.csv')

per_condition_subset_results = read_csv('../../../../intermediate/per_condition_subset_results.csv')
wt_results = read_csv('../../../../intermediate/wt_results.csv')
library(modelr)
per_condition_subset_results %>%
#per_condition_subset_results %>% filter(name %in% (modules %>% filter(module == 'HOT1') %>% pull(name))) %>%
left_join(wt_results %>%
dplyr::select(condition, log2FoldChange, name, padj) %>%
dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
dplyr::filter(wtp<0.1) %>% group_by(condition, Kinase) %>% nest() %>%
mutate(mod = map(data,~lm(log2FoldChange ~ wt_change, data = .)),
coefs = map(mod, ~coef(.)),
slope = map_dbl(coefs, 'wt_change'),
r2 = map2_dbl(mod, data, ~rsquare(.x,.y)),
var = map_dbl(data, function(x) x %>% pull(log2FoldChange) %>% var)) -> fits

growth = read_csv('~/Desktop/growth_rates.csv')

growth %>%
  select(Kinase = Strain,
         condition = Condition,
         growth = top_speed_n) %>%
  left_join(fits %>%
              select(Kinase, condition, slope)) %>%
  na.omit() %>%
  filter(condition == 'Salt') %>%
  ggplot(aes(x = slope, y = growth, label = Kinase)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggrepel::geom_text_repel()
ggsave('~/Desktop/weak_slope_to_growth_corrolation.pdf')



normalize_it = function(x) x/max(x)

# Average out the condition replicates
dat = dat %>% group_by(Condition, Strain) %>% summarize(top_speed = mean(top_speed)) %>% dplyr::ungroup()

dat = dat %>%
        group_by(Strain) %>% mutate(ts_st = normalize_it(top_speed)) %>% dplyr::ungroup() %>%
        group_by(Condition) %>% mutate(ts_cnd = normalize_it(top_speed)) %>% dplyr::ungroup()

# dat = dat %>%
#         group_by(Strain) %>% mutate(ts_st = normalize_it(top_speed)) %>% dplyr::ungroup() %>%
#         group_by(Condition) %>% mutate(ts_cnd = normalize_it(top_speed)) %>% dplyr::ungroup()

dat = dat %>% group_by(Strain) %>% mutate(ts_st_d = top_speed/mean(top_speed[Condition=='Drug']))

pdf('../../output/growth.R/top_speed_by_condition.pdf')
ggplot(dat, aes(x = reorder(Condition, top_speed, FUN=median), y = top_speed)) +
geom_boxplot() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
xlab("Condition") +
ylab("Fastest Growth Rate (Doublings/hr)") +
ggtitle("Growth Rate grouped by Condition")
dev.off()

pdf('../../output/growth.R/top_speed_by_strain.pdf')
ggplot(dat, aes(x = reorder(Strain, top_speed, FUN=median), y = top_speed)) +
geom_boxplot() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
xlab("Condition") +
ylab("Fastest Growth Rate (Doublings/hr)") +
ggtitle("Growth Rate by grouped by Strain")
dev.off()

pdf('../../output/growth.R/top_speed_facet_condition_same_axis.pdf')
ggplot(dat, aes(x = reorder(Strain, top_speed, FUN=median), y = top_speed)) +
geom_bar(stat="identity") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
facet_wrap(~ Condition, scale='free_x') +
xlab("Condition") +
ylab("Fastest Growth Rate (Doublings/hr)") +
ggtitle("Growth Rate by Condition")
dev.off()

pdf('../../output/growth.R/top_speed_facet_condition.pdf')
ggplot(dat, aes(x = reorder(Strain, top_speed, FUN=median), y = top_speed)) +
geom_bar(stat="identity") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
facet_wrap(~ Condition, scale='free') +
xlab("Condition") +
ylab("Fastest Growth Rate (Doublings/hr)") +
ggtitle("Growth Rate by Condition")
dev.off()

pdf('../../output/growth.R/top_speed_facet_strain.pdf')
ggplot(dat, aes(x = reorder(Condition, top_speed, FUN=median), y = top_speed)) +
geom_bar(stat="identity") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
facet_wrap(~ Strain) +
xlab("Condition") +
ylab("Fastest Growth Rate (Doublings/hr)") +
ggtitle("Growth Rate by Strain")
dev.off()

pdf('../../output/growth.R/top_speed_facet_strain_no_ypd.pdf')
ggplot(filter(dat,Condition!='YPD'), aes(x = reorder(Condition, top_speed, FUN=median), y = top_speed)) +
geom_bar(stat="identity") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.4)) +
facet_wrap(~ Strain) +
xlab("Condition") +
ylab("Fastest Growth Rate (Doublings/hr)") +
ggtitle("Growth Rate by Strain")
dev.off()

# pdf('../../output/growth.R/top_speed_facet_condition_norm.pdf')
# ggplot(dat, aes(x = reorder(Strain, ts_cnd, FUN=median), y = ts_cnd)) +
# geom_bar(stat="identity") +
# theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
# facet_wrap(~ Condition, scale='free') +
# xlab("Condition") +
# ylab("Fastest Growth Rate (Doublings/hr), relative to general condition performance") +
# ggtitle("Relative Growth Rate by Condition")
# dev.off()

pdf('../../output/growth.R/top_speed_facet_condition_norm.pdf')
ggplot(dat, aes(x = reorder(Strain, ts_st, FUN=median), y = ts_st)) +
geom_bar(stat="identity") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
facet_wrap(~ Condition, scale='free') +
xlab("Condition") +
ylab("Fastest Growth Rate (Doublings/hr), relative to general condition performance") +
ggtitle("Relative Growth Rate by Condition")
dev.off()

pdf('../../output/growth.R/top_speed_facet_strain_drug_norm.pdf')
ggplot(dat, aes(x = reorder(Strain, ts_st_d, FUN=median), y = ts_st_d)) +
geom_bar(stat="identity") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
facet_wrap(~ Condition, scale='free') +
xlab("Condition") +
ylab("Fastest Growth Rate (Doublings/hr), relative to drug performance") +
ggtitle("Relative Growth Rate by Condition")
dev.off()

pdf('../../output/growth.R/top_speed_facet_strain_norm.pdf')
ggplot(dat, aes(x = reorder(Condition, ts_st, FUN=median), y = ts_st)) +
geom_bar(stat="identity") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
facet_wrap(~ Strain) +
xlab("Condition") +
ylab("Fastest Growth Rate (Doublings/hr), relative to general strain performance") +
ggtitle("Growth Rate by Strain")
dev.off()
