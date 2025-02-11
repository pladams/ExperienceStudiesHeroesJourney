library(data.table)
library(tidyverse)

nHorizon <- 3
sFracType <- "CF"
nSims <- 100000

# First decrement
annual_rate <- 0.5
daily_rate <- 1-(1-annual_rate)^(1/365)

# Second decrement
annual_rate2 <- 0.25
modal_rate2 <- 1-(1-annual_rate2)^(1/12)

set.seed(0xBEEF)

source("lifetable_to_distribution.R")

death_table <- lifetable_to_distribution(rep(daily_rate,365*nHorizon))

events <- death_table[,sample(x=time,size=nSims,replace=T,prob=q_xpt)]

lapse_table <- lifetable_to_distribution(rep(modal_rate2,nHorizon*12))

events2 <- lapse_table[,sample(x=time,size=nSims,replace=T,prob=q_xpt)]

bExactTerminations <- FALSE


study_start <- as_datetime(ymd(20210101))
study_end <- study_start %m+% years(nHorizon) %m+% days(-1)
if(bExactTerminations) {
  study_end <- study_end %m+% seconds(86400-1)
}


pol_period_granularity <- 12 # months per policy period
cal_period_granularity <- 12 # months per calendar period

cal_yr_breaks <- (study_start %m+% days(-1)) %m+% 
  months( 
    cal_period_granularity*(1:(interval(study_start %m+% days(-1),study_end) %/% months(cal_period_granularity)))
  )

if(bExactTerminations) {
  cal_yr_breaks <- cal_yr_breaks %m+% seconds(86400-1)
}

data.table(
  PolID=1:100000,
  Issue_Date=study_start %m+% days(sample.int(365,100000,replace=T)),
  Prem_Mode_Months=12
) %>%
  mutate(  Death_Date=Issue_Date %m+% days(events),
           Lapse_Date=Issue_Date %m+% months(events2)) -> 
  census

if(bExactTerminations) {
  census %>%
    mutate(
      Death_Date = Death_Date %m+% seconds(86400-1) %m+% seconds(-sample.int(86400-1,nrow(.),replace=T)),
      Lapse_Date = Lapse_Date %m+% seconds(86400-1)
    ) ->
    census
}

census %>%
  mutate(Term_Date = as_datetime(ifelse(Death_Date <= Lapse_Date & Death_Date <= study_end, Death_Date,NA))) %>%
  mutate(Term_Date = as_datetime(ifelse(Lapse_Date <= Death_Date & Lapse_Date <= study_end, Lapse_Date,Term_Date))) -> 
  census


source('expand.exposures.R')

census %>%
  filter(Issue_Date <= study_end) %>%
  select(PolID,Issue_Date,Term_Date,Prem_Mode_Months) %>%
  expand_exposures(
    .exp_period_start = study_start,
    .exp_period_end = study_end,
    .cal_yr_breaks = cal_yr_breaks,
    .pol_period_granularity = pol_period_granularity,
    .cal_period_granularity = cal_period_granularity,
    .issue_date = Issue_Date,
    .term_date = Term_Date,
    .ID=PolID,
    .prem_mode_months = Prem_Mode_Months,
    .exact_terminations = FALSE
  ) ->
  exposures

census %>%
  inner_join(y=exposures,by=join_by(PolID)) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0),
         Expected_Deaths=-(exposure)*log(1-annual_rate),
         Expected_Deaths2=annual_rate*(exposure),
         Expected_Deaths3=annual_rate*(exposure + Death_Count/2),
         Expected_Deaths4=exposure*annual_rate/(1-exposure*annual_rate)) %>%
  summarize(Death_Count=sum(Death_Count),
            Exposure=sum(exposure),
            Expected_Deaths=sum(Expected_Deaths),
            Expected_Deaths2=sum(Expected_Deaths2),
            Expected_Deaths3=sum(Expected_Deaths3),
            Expected_Deaths4=sum(Expected_Deaths4),
            Hazard=Death_Count/sum(exposure),
            Rate=1-exp(-Hazard),
            Rate2=Death_Count/(Exposure+Death_Count/2),
            A_E_1=Death_Count/Expected_Deaths,
            A_E_2=Death_Count/Expected_Deaths2,
            A_E_3=Death_Count/Expected_Deaths3,
            A_E_4=Death_Count/Expected_Deaths4,
            n()
  )



death_table_a <- data.table(
  ID=1:(nHorizon*365),
  Rate=1-(1-0.4)^(1/365)
)

death_table_a[,Rate:=Rate*1.1^((ID-1)/365)]

death_table_a <- lifetable_to_distribution(death_table_a[,Rate])
events_a <- death_table_a[,sample(x=time,size=nSims,replace=T,prob=q_xpt)]


copy(census) %>%
  mutate(  Death_Date=Issue_Date %m+% days(events_a)) -> 
  census2

if(bExactTerminations) {
  census2 %>%
    mutate(
      Death_Date = Death_Date %m+% seconds(86400-1) %m+% seconds(-sample.int(86400-1,nrow(.),replace=T)),
      Lapse_Date = Lapse_Date %m+% seconds(86400-1)
    ) ->
    census2
}

census2 %>%
  mutate(Term_Date = as_datetime(ifelse(Death_Date <= Lapse_Date & Death_Date <= study_end, Death_Date,NA))) %>%
  mutate(Term_Date = as_datetime(ifelse(Lapse_Date <= Death_Date & Lapse_Date <= study_end, Lapse_Date,Term_Date))) -> 
  census2


census2 %>%
  filter(Issue_Date <= study_end) %>%
  select(PolID,Issue_Date,Term_Date,Prem_Mode_Months) %>%
  expand_exposures(
    .exp_period_start = study_start,
    .exp_period_end = study_end,
    .cal_yr_breaks = cal_yr_breaks,
    .pol_period_granularity = pol_period_granularity,
    .cal_period_granularity = cal_period_granularity,
    .issue_date = Issue_Date,
    .term_date = Term_Date,
    .ID=PolID,
    .prem_mode_months = Prem_Mode_Months,
    .exact_terminations = FALSE
  ) ->
  exposures2


census2 %>%
  inner_join(y=exposures2,by=join_by(PolID)) %>%
  inner_join(
    y=death_table_a %>%
      filter(time %in% c(365*(1:3))) %>%
      mutate(
        pol_duration=time/365,
        qx=1-tpx/shift(tpx,fill=1)
      ) %>%
      select(
        pol_duration,
        qx
      ),
    by=join_by(pol_duration)
  ) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0),
         Expected_Deaths=-(exposure)*log(1-qx),
         Expected_Deaths2=qx*(exposure),
         Expected_Deaths3=qx*(exposure + Death_Count/2),
         Expected_Deaths4=exposure*qx/(1-exposure*qx)) %>%
  group_by(pol_duration) %>%
  summarize(Death_Count=sum(Death_Count),
            Exposure=sum(exposure),
            Hazard=sum(Death_Count)/sum(exposure),
            Rate=1-exp(-Hazard),
            Rate2=sum(Death_Count)/sum((Exposure+Death_Count/2)),
            A_E_1=sum(Death_Count)/sum(Expected_Deaths),
            A_E_2=sum(Death_Count)/sum(Expected_Deaths2),
            A_E_3=sum(Death_Count)/sum(Expected_Deaths3),
            A_E_4=sum(Death_Count)/sum(Expected_Deaths4)
  )

death_table_a %>%
  filter(time %in% c(round(365*(1:6)/2,0))) %>%
  mutate(
    pol_duration=time/365,
    qx=1-tpx/shift(tpx,fill=1)
  ) %>%
  mutate(
    qx=1-(1-qx)^2
  ) %>%
  select(
    pol_duration,
    qx
  )


copy(census) %>%
  mutate(  Death_Date=Issue_Date %m+% days(ceiling(events/365)*365)) -> 
  census3

if(bExactTerminations) {
  census3 %>%
    mutate(
      Death_Date = Death_Date %m+% seconds(86400-1) %m+% seconds(-sample.int(86400-1,nrow(.),replace=T)),
      Lapse_Date = Lapse_Date %m+% seconds(86400-1)
    ) ->
    census3
}

census3 %>%
  mutate(Term_Date = as_datetime(ifelse(Death_Date <= Lapse_Date & Death_Date <= study_end, Death_Date,NA))) %>%
  mutate(Term_Date = as_datetime(ifelse(Lapse_Date <= Death_Date & Lapse_Date <= study_end, Lapse_Date,Term_Date))) -> 
  census3


census3 %>%
  filter(Issue_Date <= study_end) %>%
  select(PolID,Issue_Date,Term_Date,Prem_Mode_Months) %>%
  expand_exposures(
    .exp_period_start = study_start,
    .exp_period_end = study_end,
    .cal_yr_breaks = cal_yr_breaks,
    .pol_period_granularity = pol_period_granularity,
    .cal_period_granularity = cal_period_granularity,
    .issue_date = Issue_Date,
    .term_date = Term_Date,
    .ID=PolID,
    .prem_mode_months = Prem_Mode_Months,
    .exact_terminations = FALSE
  ) ->
  exposures3

census3 %>%
  mutate(  Death_Date=Issue_Date %m+% days(events)) -> 
  census3


census3 %>%
  inner_join(y=exposures3,by=join_by(PolID)) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0),
         Expected_Deaths=-(exposure)*log(1-annual_rate),
         Expected_Deaths2=annual_rate*(exposure),
         Expected_Deaths3=annual_rate*(exposure + Death_Count/2),
         Expected_Deaths4=exposure*annual_rate/(1-exposure*annual_rate)) %>%
  group_by(Year=year(exp_period_end)) %>%
  summarize(Death_Count=sum(Death_Count),
            Exposure=sum(exposure),
            Hazard=sum(Death_Count)/sum(exposure),
            Rate=1-exp(-Hazard),
            Rate2=sum(Death_Count)/sum((Exposure+Death_Count/2)),
            A_E_1=sum(Death_Count)/sum(Expected_Deaths),
            A_E_2=sum(Death_Count)/sum(Expected_Deaths2),
            A_E_3=sum(Death_Count)/sum(Expected_Deaths3),
            A_E_4=sum(Death_Count)/sum(Expected_Deaths4),
            n()
  )

  


copy(census2) %>%
  mutate(  Death_Date=Issue_Date %m+% days(ceiling(events_a/365)*365)) -> 
  census4

if(bExactTerminations) {
  census4 %>%
    mutate(
      Death_Date = Death_Date %m+% seconds(86400-1) %m+% seconds(-sample.int(86400-1,nrow(.),replace=T)),
      Lapse_Date = Lapse_Date %m+% seconds(86400-1)
    ) ->
    census4
}

census4 %>%
  mutate(Term_Date = as_datetime(ifelse(Death_Date <= Lapse_Date & Death_Date <= study_end, Death_Date,NA))) %>%
  mutate(Term_Date = as_datetime(ifelse(Lapse_Date <= Death_Date & Lapse_Date <= study_end, Lapse_Date,Term_Date))) -> 
  census4


census4 %>%
  filter(Issue_Date <= study_end) %>%
  select(PolID,Issue_Date,Term_Date,Prem_Mode_Months) %>%
  expand_exposures(
    .exp_period_start = study_start,
    .exp_period_end = study_end,
    .cal_yr_breaks = cal_yr_breaks,
    .pol_period_granularity = pol_period_granularity,
    .cal_period_granularity = cal_period_granularity,
    .issue_date = Issue_Date,
    .term_date = Term_Date,
    .ID=PolID,
    .prem_mode_months = Prem_Mode_Months,
    .exact_terminations = FALSE
  ) ->
  exposures4

census4 %>%
  mutate(  Death_Date=Issue_Date %m+% days(events_a)) -> 
  census4


census4 %>%
  inner_join(y=exposures4,by=join_by(PolID)) %>%
  inner_join(
    y=death_table_a %>%
      filter(time %in% c(365*(1:3))) %>%
      mutate(
        pol_duration=time/365,
        qx=1-tpx/shift(tpx,fill=1)
      ) %>%
      select(
        pol_duration,
        qx
      ),
    by=join_by(pol_duration)
  ) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0),
         Expected_Deaths=-(exposure)*log(1-qx),
         Expected_Deaths2=qx*(exposure),
         Expected_Deaths3=qx*(exposure + Death_Count/2),
         Expected_Deaths4=exposure*qx/(1-exposure*qx)) %>%
  group_by(pol_duration) %>%
  summarize(Death_Count=sum(Death_Count),
            Exposure=sum(exposure),
            Hazard=sum(Death_Count)/sum(exposure),
            Rate=1-exp(-Hazard),
            Rate2=sum(Death_Count)/sum((Exposure+Death_Count/2)),
            A_E_1=sum(Death_Count)/sum(Expected_Deaths),
            A_E_2=sum(Death_Count)/sum(Expected_Deaths2),
            A_E_3=sum(Death_Count)/sum(Expected_Deaths3),
            A_E_4=sum(Death_Count)/sum(Expected_Deaths4)
  )


census4 %>%
  inner_join(y=exposures4,by=join_by(PolID)) %>%
  inner_join(
    y=death_table_a %>%
      filter(time %in% c(365*(1:3))) %>%
      mutate(
        pol_duration=time/365,
        qx=1-tpx/shift(tpx,fill=1)
      ) %>%
      select(
        pol_duration,
        qx
      ),
    by=join_by(pol_duration)
  ) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0),
         Expected_Deaths=-(exposure)*log(1-qx),
         Expected_Deaths2=qx*(exposure),
         Expected_Deaths3=qx*(exposure + Death_Count/2),
         Expected_Deaths4=exposure*qx/(1-exposure*qx)) %>%
  group_by(Year=year(exp_period_end)) %>%
  summarize(Death_Count=sum(Death_Count),
            Exposure=sum(exposure),
            Hazard=sum(Death_Count)/sum(exposure),
            Rate=1-exp(-Hazard),
            Rate2=sum(Death_Count)/sum((Exposure+Death_Count/2)),
            A_E_1=sum(Death_Count)/sum(Expected_Deaths),
            A_E_2=sum(Death_Count)/sum(Expected_Deaths2),
            A_E_3=sum(Death_Count)/sum(Expected_Deaths3),
            A_E_4=sum(Death_Count)/sum(Expected_Deaths4)
  )

#GLM
census %>%
  inner_join(y=exposures,by=join_by(PolID)) %>%
  filter(exposure > 0) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0)) %>%
  glm(
    Death_Count ~ offset(log(exposure)),
    data=.,
    family=poisson
  )

census2 %>%
  inner_join(y=exposures2,by=join_by(PolID)) %>%
  filter(exposure > 0) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0)) %>%
  glm(
    Death_Count ~ factor(pol_duration) + offset(log(exposure)),
    data=.,
    family=poisson
  )


death_table_a %>%
  filter(time %in% c(365*(1:3))) %>%
  mutate(
    pol_duration=time/365,
    qx=1-tpx/shift(tpx,fill=1)
  ) %>%
  select(
    pol_duration,
    qx
  ) %>%
  mutate(
    Hazard=-log(1-qx)
  )


census3 %>%
  inner_join(y=exposures3,by=join_by(PolID)) %>%
  filter(exposure > 0) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0)) %>%
  glm(
    Death_Count ~ offset(log(exposure)),
    data=.,
    family=poisson
  )


census4 %>%
  inner_join(y=exposures4,by=join_by(PolID)) %>%
  filter(exposure > 0) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0)) %>%
  glm(
    Death_Count ~ factor(pol_duration) + offset(log(exposure)),
    data=.,
    family=poisson
  )


census4 %>%
  inner_join(y=exposures4,by=join_by(PolID)) %>%
  filter(exposure > 0) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0)) %>%
  glm(
    Death_Count ~ offset(log(exposure)),
    data=.,
    family=poisson
  )
