# create shiny.RData
load('../../../intermediate/images/paper_data.RData')

genes %>%
  select(name, data) %>%
  unnest() %>%
  select(name, Strain_Code, Condition, Differential_Expression) %>%
  left_join(modules) %>%
  group_by(module, Strain_Code, Condition) %>%
  summarize(avg_d_exp = mean(Differential_Expression)) %>%
  ungroup() %>% rename(Kinase = Strain_Code) -> module_expression

tfs = module_tfs

save(list = c(
"genes",
"MEs",
"MEsCor",
"module_raw_residuals",
"module_std_residuals",
"module_expression",
"modules",
"raw_residuals",
"std_residuals",
"tfs"),
file = 'shiny.RData')

