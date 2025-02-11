


# Validating sampler
test.samples <- replicate(
  n=1000000,
  expr=sample_death_duration("F",55,2045,1)
)

saveRDS(test.samples,"test.samples.rds")

test.samples %>%
  data.table(Duration=.) %>%
  group_by(Duration) %>%
  summarize(Deaths=n()) %>%
  ungroup() %>%
  right_join(y=data.table(Duration=1:64)) %>%
  mutate(Deaths=ifelse(is.na(Deaths),0,Deaths)) %>%
  arrange(-Duration) %>%
  mutate(Exposures=cumsum(Deaths)) %>%
  arrange(Duration) %>%
  mutate(qx_sim=Deaths/Exposures) %>%
  inner_join(
    y=vbt2015[Sex=="F" & Issue_Age==55]
  ) %>%
  mutate(ratio=qx_sim/qx_su,
         se2=2/sqrt(Deaths),
         pearson_resid=(Deaths-Exposures*qx_su)/sqrt(Exposures*qx_su)) %>%
  data.table() %>%
  ggplot(aes(x=Duration,y=pearson_resid)) +
  geom_line() +
  geom_hline(yintercept=c(-2,2))
