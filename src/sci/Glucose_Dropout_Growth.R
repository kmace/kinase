source('../utils/growth.R')
# Build Meta
row_meta = read.csv('../../input/meta/Glucose_Dropout_rows.csv',header=T)
col_meta = read.csv('../../input/meta/Glucose_Dropout_cols.csv',header=T)
meta = merge(row_meta, col_meta)
meta = meta %>% mutate(Well = paste(Row, Column,sep=''),
                       Repeat = if_else(Row %in% c('A','B','C','D'),1,2))
meta$Condition = factor(meta$Condition, c("0%", "0.02%", "0.2%", "2%"))

# Load od data
od_left = process_plate('../../input/growth/Glucose_Dropout/Left_Shark/KGM_PostStress_20170316_134258.xlsx',meta)
od_right = process_plate('../../input/growth/Glucose_Dropout/Right_Shark/KGM_PostStress_20170316_134439.xlsx',meta)

# This is Pre-Stress data
# od_left = process_plate('../../input/growth/Glucose_Dropout/Left_Shark/KGM_PreStress_20170315_142703.xlsx',meta)
# od_right = process_plate('../../input/growth/Glucose_Dropout/Right_Shark/KGM_PreStress_20170315_142613.xlsx',meta)

left_sum = SummarizeGrowthByPlate(od_left, bg_correct = 'blank', plot_fit=TRUE, plot_file='left.pdf') %>% rename(Well = sample)
right_sum = SummarizeGrowthByPlate(od_right, bg_correct = 'blank',  plot_fit=TRUE, plot_file='right.pdf') %>% rename(Well = sample)
all_sum = join_plates(left_sum,right_sum,meta)

left_time = get_temporal(od_left)
right_time = get_temporal(od_right)
all_time = join_plates(left_time,right_time,meta)

all = inner_join(all_time,all_sum)
all = all %>% filter(Strain != 'Blank') %>% mutate(Predicted_OD = NAtT(k, n0, r, Time))
all$Repeat = factor(all$Repeat)

all %>% ggplot(aes(Drug, r)) + geom_boxplot() + geom_point() + facet_grid(Condition~Strain)

all %>% ggplot(aes(x = Condition, y=r)) + geom_point() + geom_boxplot()+ facet_grid(Drug~Strain)

# Make sense, population fitness
all %>% ggplot(aes(Condition, auc_l)) + geom_boxplot() + geom_point() + facet_grid(Drug~Strain)

all %>% group_by(Condition,Strain,Drug) %>% summarize(r = mean(r), rsd = sd(r)) %>% ggplot(aes(x = Condition, y=r, color=Drug)) + geom_point() + geom_line() + facet_grid(~Strain) + ylim(0,1)








all %>% filter(Growth_Rate < 1.25) %>% na.omit() %>% ggplot(aes(x = Time, y = Growth_Rate, color=Drug, shape=Repeat)) + geom_point() + geom_line() + facet_grid(Condition~Strain)

all %>% filter() %>% na.omit() %>% ggplot(aes(x = Time, y = OD, color=Drug, shape=Repeat)) + geom_point() + geom_line() + facet_grid(Condition~Strain)

all %>% na.omit() %>% ggplot(aes(x = Time, y = OD, color=Drug, shape=Repeat)) + geom_point() + geom_line() + facet_grid(Condition~Strain)

all %>% na.omit() %>% ggplot(aes(x = Time, y = OD, color=Condition, group=paste0(Repeat,Drug,Condition), shape = Drug)) + geom_point() + geom_line() + facet_wrap(~Strain)

all %>% na.omit() %>% ggplot(aes(x = Time, y = Growth_Rate, color=Condition, group=paste0(Repeat,Drug,Condition), shape = Drug)) + geom_point() + geom_line() + facet_wrap(~Strain)
