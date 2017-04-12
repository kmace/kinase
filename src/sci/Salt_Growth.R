library(dplyr)
library(tidyr)
library(ggplot2)
meta = read.csv('../../input/meta/Salt_plate_meta.csv',header=T)
od = read.csv('../../input/growth/salt.csv', header=T)
colnames(od)[1:3] = c("Cycle", 'Time', 'Temp')


rate = 3*apply(log2(apply(od[,-c(1,2,3)],2,smooth)),2,diff)
rate = rbind(rate, rep(NA, dim(rate)[2]))
rate = cbind(od[,1:3], rate)
meta = meta %>% mutate(Well = paste(Row, Column,sep=''))
meta$Repeat = factor(meta$Repeat)
timestep = 20 # min
#od = od[1:45,]

od = od %>% gather(Well, OD, -Temp, -Time, -Cycle)
rate = rate %>% gather(Well, Growth_Rate, -Temp, -Time, -Cycle)
all = left_join(od, rate)
all = left_join(all, meta)

all = all %>% filter(Column!=3 & Column!=5) # Messed up these wells


all %>% na.omit() %>% ggplot(aes(x = Time, y = OD, color=Condition, shape=Repeat)) + geom_point() + geom_line() + facet_wrap(~Strain)

all %>% na.omit() %>% ggplot(aes(x = Time, y = Growth_Rate, color=Condition, shape=Repeat)) + geom_point() + geom_line() + facet_wrap(~Strain)

p = all %>% group_by(Strain, Condition, Time) %>% summarise(OD = mean(OD)) %>%ggplot(aes(x = Time, y = OD, color = Condition)) + geom_line() + facet_wrap(~Strain)
p
ggsave('salt_growth_od.pdf')

fastest_wt = all %>%
             filter(grepl("WT", Strain)) %>%
             group_by(Condition, Cycle) %>%
             summarize(rate = mean(Growth_Rate)) %>%
             ungroup() %>%
             group_by(Condition) %>%
             slice(which.max(rate))

only_best = semi_join(all, fastest_wt)

only_best %>%
group_by(Strain, Condition) %>%
summarise(OD = mean(OD),
          Growth_Rate = mean(Growth_Rate)) %>%
ungroup() %>%
ggplot(aes(x=Condition, y = Growth_Rate)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Strain) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


only_best %>%
group_by(Strain, Condition) %>%
summarise(Growth_Rate = mean(Growth_Rate)) %>%
ungroup() %>% group_by(Strain) %>%
mutate(GRN = Growth_Rate / Growth_Rate[Condition=='Salt']) %>%
ggplot(aes(x = Condition, y = (GRN))) + geom_bar(stat='identity') +
facet_wrap(~Strain) +
theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calculate Parameters:
# Fastest Growth Rate
ave_top_speed = apply(rate,2,function(x) mean(sort(x,decreasing=TRUE)[1:3]))
dat = cbind(meta, ave_top_speed[match(meta$well_name, names(ave_top_speed))])
colnames(dat)[dim(dat)[2]] = 'top_speed'
# Time to Fastest Growth Rate
time_to_top_speed = unlist(apply(rate,2,function(x) which(x==max(x))))
dat = cbind(dat, time_to_top_speed[match(dat$well_name, names(time_to_top_speed))])
colnames(dat)[dim(dat)[2]] = 'time_to_top_speed'

# Average out the condition replicates
dat = dat %>% group_by(Condition, Strain) %>%
summarize(top_speed = mean(top_speed),
          time_to_top_speed = mean(time_to_top_speed)) %>%
          dplyr::ungroup()


normalize_it = function(x) x/max(x)

# dat = dat %>%
#   group_by(Strain) %>% mutate(ts_st = normalize_it(top_speed)) %>% dplyr::ungroup() %>%
#   group_by(Condition) %>% mutate(ts_cnd = normalize_it(top_speed)) %>% dplyr::ungroup()
#
# # dat = dat %>%
#         group_by(Strain) %>% mutate(ts_st = normalize_it(top_speed)) %>% dplyr::ungroup() %>%
#         group_by(Condition) %>% mutate(ts_cnd = normalize_it(top_speed)) %>% dplyr::ungroup()

# Strain Specific Measurements
dat = dat %>% group_by(Strain) %>%
mutate(top_speed_rel_drug = top_speed/mean(top_speed[Condition=='Drug']),
       top_speed_rel_salt = top_speed/mean(top_speed[Condition=='Salt']),
       top_speed_rel_ypd = top_speed/mean(top_speed[Condition=='YPD']),
       time_to_top_speed_rel_drug = time_to_top_speed/mean(time_to_top_speed[Condition=='Drug']),
       time_to_top_speed_rel_salt = time_to_top_speed/mean(time_to_top_speed[Condition=='Salt']),
       time_to_top_speed_rel_ypd = time_to_top_speed/mean(time_to_top_speed[Condition=='YPD']))


dir.create('../../output/Salt_Growth.R')

pdf('../../output/Salt_Growth.R/growth_by_strain.pdf')
for(s in sort(as.character(unique(meta$Strain)))){
  strain = s;
  m = filter(meta, Strain==strain);
  matplot((0:(dim(rate)[1]-1) * timestep)/60,
          apply(rate[,m$well_name],2,smooth),
          type='l',
          col=brewer.pal(8, "Set2")[factor(m$Condition)],
          lty='solid',
          ylab='Doublings per hour',
          ylim=c(-0.1,1),
          lwd=1.5,
          xlab='Time (h)',
          main = strain);
  legend('topright',legend=m$Condition,
         col=brewer.pal(8, "Set2")[factor(m$Condition)],
         lty='solid',
         lwd=6)
}
dev.off()


pdf('../../output/Salt_Growth.R/growths_by_condition.pdf')
for(s in as.character(unique(meta$Condition))){
  condition = s;
  m = filter(meta, Condition==condition);
  matplot((0:(dim(rate)[1]-1) * timestep)/60,
          apply(rate[,m$well_name],2,smooth),
          type='l',
          col=brewer.pal(11, "Set3")[as.numeric(factor(m$Strain)) %% 11 + 1],
          lty=c('solid','dashed', 'dotted')[as.numeric(factor(m$Strain)) %% 3 + 1],
          lwd=1.5,
          ylab='Doublings per hour',
          ylim=c(-0.1,1),
          xlab='Time (h)',
          main = condition);
  legend('topleft',legend=m$Strain,
         col=brewer.pal(11, "Set3")[as.numeric(factor(m$Strain)) %% 11 + 1],
         lty=c('solid','dashed', 'dotted')[as.numeric(factor(m$Strain)) %% 3 + 1],
         lwd=1.5,
         cex=0.5)
}
dev.off()




# pdf('../../output/Salt_Growth.R/top_speed_by_condition.pdf')
# ggplot(dat, aes(x = reorder(Condition, top_speed, FUN=median), y = top_speed)) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
#   xlab("Condition") +
#   ylab("Fastest Growth Rate (Doublings/hr)") +
#   ggtitle("Growth Rate grouped by Condition")
# dev.off()

# pdf('../../output/Salt_Growth.R/top_speed_by_strain.pdf')
# ggplot(dat, aes(x = reorder(Strain, top_speed, FUN=median), y = top_speed)) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
#   xlab("Condition") +
#   ylab("Fastest Growth Rate (Doublings/hr)") +
#   ggtitle("Growth Rate by grouped by Strain")
# dev.off()

# pdf('../../output/Salt_Growth.R/top_speed_facet_condition.pdf')
# ggplot(dat, aes(x = reorder(Strain, top_speed, FUN=median), y = top_speed)) +
#   geom_bar(stat="identity") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
#   facet_wrap(~ Condition, scale='free') +
#   xlab("Condition") +
#   ylab("Fastest Growth Rate (Doublings/hr)") +
#   ggtitle("Growth Rate by Condition")
# dev.off()

# Just make all the Y axis the same for comaprision across conditions

pdf('../../output/Salt_Growth.R/top_speed_facet_condition.pdf')
ggplot(dat, aes(x = Strain, y = top_speed)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
  facet_wrap(~ Condition, scale='free_x') +
  xlab("Condition") +
  ylab("Fastest Growth Rate (Doublings/hr)") +
  ggtitle("Growth Rate by Condition")
dev.off()

# Across Strains
pdf('../../output/Salt_Growth.R/top_speed_facet_strain.pdf')
ggplot(dat, aes(x = Condition, y = top_speed)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
  facet_wrap(~ Strain) +
  xlab("Condition") +
  ylab("Fastest Growth Rate (Doublings/hr)") +
  ggtitle("Growth Rate by Strain")
dev.off()

pdf('../../output/Salt_Growth.R/time_to_top_speed_facet_condition.pdf')
ggplot(dat, aes(x = Strain, y = time_to_top_speed)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
  facet_wrap(~ Condition, scale='free_x') +
  xlab("Condition") +
  ylab("Fastest Growth Rate (Doublings/hr)") +
  ggtitle("Growth Rate by Condition")
dev.off()

# Across Strains
pdf('../../output/Salt_Growth.R/time_to_top_speed_facet_strain.pdf')
ggplot(dat, aes(x = Condition, y = time_to_top_speed)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
  facet_wrap(~ Strain) +
  xlab("Condition") +
  ylab("Fastest Growth Rate (Doublings/hr)") +
  ggtitle("Growth Rate by Strain")
dev.off()

# pdf('../../output/Salt_Growth.R/top_speed_facet_strain_no_ypd.pdf')
# ggplot(filter(dat,Condition!='YPD'), aes(x = reorder(Condition, top_speed, FUN=median), y = top_speed)) +
#   geom_bar(stat="identity") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.4)) +
#   facet_wrap(~ Strain) +
#   xlab("Condition") +
#   ylab("Fastest Growth Rate (Doublings/hr)") +
#   ggtitle("Growth Rate by Strain")
# dev.off()

# pdf('../../output/Salt_Growth.R/top_speed_facet_condition_norm.pdf')
# ggplot(dat, aes(x = reorder(Strain, ts_cnd, FUN=median), y = ts_cnd)) +
# geom_bar(stat="identity") +
# theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
# facet_wrap(~ Condition, scale='free') +
# xlab("Condition") +
# ylab("Fastest Growth Rate (Doublings/hr), relative to general condition performance") +
# ggtitle("Relative Growth Rate by Condition")
# dev.off()

# pdf('../../output/Salt_Growth.R/top_speed_facet_condition_norm.pdf')
# ggplot(dat, aes(x = reorder(Strain, ts_st, FUN=median), y = ts_st)) +
#   geom_bar(stat="identity") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
#   facet_wrap(~ Condition, scale='free') +
#   xlab("Condition") +
#   ylab("Fastest Growth Rate (Doublings/hr), relative to general condition performance") +
#   ggtitle("Relative Growth Rate by Condition")
# dev.off()

# Max Growht

pdf('../../output/Salt_Growth.R/top_speed_facet_stress_drug_norm.pdf')
ggplot(dat, aes(x = Strain, y = top_speed_rel_drug)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
  facet_wrap(~ Condition, scale='free') +
  xlab("Condition") +
  ylab("Fastest Growth Rate (Doublings/hr), relative to drug performance") +
  ggtitle("Relative Growth Rate by Condition")
dev.off()

pdf('../../output/Salt_Growth.R/top_speed_facet_stress_stress_norm.pdf')
ggplot(dat, aes(x = Strain, y = top_speed_rel_salt)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
  facet_wrap(~ Condition, scale='free') +
  xlab("Condition") +
  ylab("Fastest Growth Rate (Doublings/hr), relative to stress performance") +
  ggtitle("Relative Growth Rate by Condition")
dev.off()

pdf('../../output/Salt_Growth.R/top_speed_facet_stress_ypd_norm.pdf')
ggplot(dat, aes(x = Strain, y = top_speed_rel_ypd)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
  facet_wrap(~ Condition, scale='free') +
  xlab("Condition") +
  ylab("Fastest Growth Rate (Doublings/hr), relative to stress performance") +
  ggtitle("Relative Growth Rate by Condition")
dev.off()

# Time to max growth

pdf('../../output/Salt_Growth.R/time_to_top_speed_facet_stress_drug_norm.pdf')
ggplot(dat, aes(x = Strain, y = time_to_top_speed_rel_drug)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
  facet_wrap(~ Condition, scale='free') +
  xlab("Condition") +
  ylab("Fastest Growth Rate (Doublings/hr), relative to drug performance") +
  ggtitle("Relative Growth Rate by Condition")
dev.off()

pdf('../../output/Salt_Growth.R/time_to_top_speed_facet_stress_stress_norm.pdf')
ggplot(dat, aes(x = Strain, y = time_to_top_speed_rel_salt)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
  facet_wrap(~ Condition, scale='free') +
  xlab("Condition") +
  ylab("Fastest Growth Rate (Doublings/hr), relative to stress performance") +
  ggtitle("Relative Growth Rate by Condition")
dev.off()

pdf('../../output/Salt_Growth.R/time_to_top_speed_facet_stress_ypd_norm.pdf')
ggplot(dat, aes(x = Strain, y = time_to_top_speed_rel_ypd)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 4)) +
  facet_wrap(~ Condition, scale='free') +
  xlab("Condition") +
  ylab("Fastest Growth Rate (Doublings/hr), relative to stress performance") +
  ggtitle("Relative Growth Rate by Condition")
dev.off()

# pdf('../../output/Salt_Growth.R/top_speed_facet_strain_norm.pdf')
# ggplot(dat, aes(x = reorder(Condition, ts_st, FUN=median), y = ts_st)) +
#   geom_bar(stat="identity") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
#   facet_wrap(~ Strain) +
#   xlab("Condition") +
#   ylab("Fastest Growth Rate (Doublings/hr), relative to general strain performance") +
#   ggtitle("Growth Rate by Strain")
# dev.off()
