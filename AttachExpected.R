library(data.table)
library(tidyverse)
library(visreg)


# Develop pricing assumptions, based on issue years 2041-2045, experiene year 2041-2045

# Generate (simple) mortality assumption for pricing
policy_pop %>%
  filter(UW_Decision == "STD" & year(Issue_Date) <= 2046) %>%
  inner_join(y=policy_exposures,
             by="PolID") %>%
  filter( year(exp_period_end) <= 2046) %>%
  inner_join(
    y=vbt2015,
    by=join_by(
      Sex==Sex,
      Issue_Age==Issue_Age,
      pol_duration==Duration
    )
  ) %>%
  mutate(ExpectedClaims_2015VBT_Count=exposure*-log(1-qx_su),
         ExpectedClaims_2015VBT_Amount=ExpectedClaims_2015VBT_Count*Face_Amount,
         pol_duration_mort_fit=factor(pol_duration)) %>%
  group_by(Sex,PrefClass=factor(PrefClass),Issue_Age,pol_duration_mort_fit,Face_Group) %>%
  summarize(
    Death_Count=sum(ifelse(Term_Date <= exp_period_end & Term_Date >= exp_period_start & Status == "Death",1,0)),
    ExpectedClaims_2015VBT_Count=sum(ExpectedClaims_2015VBT_Count),
    Death_Amount=sum(ifelse(Term_Date <= exp_period_end & Term_Date >= exp_period_start & Status == "Death",Face_Amount,0)),
    ExpectedClaims_2015VBT_Amount=sum(ExpectedClaims_2015VBT_Amount)
  ) %>%
  data.table() ->
  pricing.mortality

glm(Death_Count ~ offset(log(ExpectedClaims_2015VBT_Count))+(Sex + PrefClass+Issue_Age + pol_duration_mort_fit + Face_Group)^2,
    family=poisson,
    data=pricing.mortality
    ) %>%
  MASS::stepAIC(
    object=.,
    test="Chisq"
  ) -> 
  mod.glm 

mod.glm %>%
  summary

mod.glm %>%
  update(object=.,
         ~.-pol_duration_mort_fit:Sex-Face_Group) ->
  mod.glm

mod.glm %>%
  summary

glm(Death_Count ~ offset(log(ExpectedClaims_2015VBT_Count))+(Sex + PrefClass+Issue_Age + I(pol_duration_mort_fit==1) + Face_Group)^2,
    family=poisson,
    data=pricing.mortality
) %>%
  MASS::stepAIC(
    object=.,
    test="Chisq"
  ) -> 
  mod.glm 

mod.glm %>%
  summary

mod.glm %>%
  update(object=.,
         ~.-I(pol_duration_mort_fit==1):Face_Group - Face_Group) ->
  mod.glm

mod.glm %>%
  summary

# The modeling has trouble with finding the subtle duration trend and the gender differential

# We ignore the credibility issues of this approach

# Lapse Assumption
policy_pop %>%
  filter(UW_Decision == "STD" & year(Issue_Date) <= 2046) %>%
  inner_join(y=policy_exposures,
             by="PolID") %>%
  filter( year(exp_period_end) <= 2046) %>%
  group_by(Sex,PrefClass=factor(PrefClass),Issue_Age,pol_duration_lapse_fit=factor(pol_duration),Face_Group) %>%
  summarize(
    Lapse_Count=sum(ifelse(Term_Date <= exp_period_end & Term_Date >= exp_period_start & Status == "Lapsed",1,0)),
    Exposure_Count=sum(exposure),
    Lapse_Amount=sum(ifelse(Term_Date <= exp_period_end & Term_Date >= exp_period_start & Status == "Lapsed",Face_Amount,0)),
    Exposure_Amount=sum(exposure*Face_Amount)
  ) %>%
  data.table()  ->
  pricing.lapse


glm(Lapse_Count ~ offset(log(Exposure_Count))+(Sex + PrefClass+Issue_Age + pol_duration_lapse_fit + Face_Group)^2,
    family=poisson,
    data=pricing.lapse
) %>%
  MASS::stepAIC(
    object=.,
    test="Chisq"
  ) -> 
  mod.glm.lapse

mod.glm.lapse %>%
  summary

mod.glm.lapse %>%
  update(object=.,
         ~.-Sex:pol_duration_lapse_fit) ->
  mod.glm.lapse

mod.glm.lapse %>%
  summary

policy_pop %>%
  filter(UW_Decision != "DEC") %>%
  inner_join(y=policy_exposures,
             by="PolID") %>%
  inner_join(
    y=vbt2015,
    by=join_by(
      Sex==Sex,
      Issue_Age==Issue_Age,
      pol_duration==Duration
    )
  ) %>%
  mutate(PrefClass=factor(PrefClass),
         Exposure_Count=exposure,
         pol_duration_lapse_fit=factor(pmin(6,pol_duration)),
         pol_duration_mort_fit=pol_duration_lapse_fit,
         ExpectedClaims_2015VBT_Count=exposure*-log(1-qx_su),
         ExpectedClaims_2015VBT_Amount=ExpectedClaims_2015VBT_Count*Face_Amount) %>%
  mutate(
    Predicted_Pricing_Mortality_Count=predict(mod.glm,newdata=.,type="response"),
    Predicted_Pricing_Lapse_Count=predict(mod.glm.lapse,newdata=.,type="response"),
    Predicted_Pricing_Mortality_Amount=Predicted_Pricing_Mortality_Count*Face_Amount,
    Predicted_Pricing_Lapse_Amount=Predicted_Pricing_Lapse_Count*Face_Amount
  ) %>%
  select(ID,
         Predicted_Pricing_Mortality_Count,
         Predicted_Pricing_Lapse_Count,
         Predicted_Pricing_Mortality_Amount,
         Predicted_Pricing_Lapse_Amount) ->
  expected_pricing

saveRDS(expected_pricing,"expected_pricing.rds")
