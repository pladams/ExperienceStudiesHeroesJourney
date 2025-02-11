library(data.table)
library(tidyverse)

policy_pop %>%
  filter(UW_Decision != "DEC") %>%
  inner_join(y=stage2,
             by="PolID") %>%
  inner_join(
    y=vbt2015,
    by=join_by(
      Sex==Sex,
      Issue_Age==Issue_Age,
      pol_duration==Duration
    )
  ) %>%
  select(-qx_ult) %>%
  mutate(ExpectedClaimsHaz_Count=exposure*-log(1-qx_su),
         ExpectedClaimsHaz_Amount=ExpectedClaimsHaz_Count*Face_Amount,
         ExpectedClaimsQx1_Count=qx_su*exposure,
         ExpectedClaimsQx2_Count=exposure*qx_su/(1-(1-exposure)*qx_su)) %>%
  select(ID,ExpectedClaimsHaz_Count,ExpectedClaimsHaz_Amount,
         ExpectedClaimsQx1_Count,ExpectedClaimsQx2_Count) ->
  expected_2015vbt

saveRDS(expected_2015vbt,"expected_2015vbt.rds")
