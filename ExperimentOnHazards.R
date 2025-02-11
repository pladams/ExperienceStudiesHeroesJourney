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
# Sample deaths
event_probs <- data.table(
  time=seq(1,nHorizon+1),
  qx=c(rep(annual_rate,nHorizon),1)
)
event_probs[,
            tpx:=cumprod(1-qx)
]
event_probs[,
            q_xpt:=qx*shift(tpx,fill=1)
]

events <- event_probs[,sample(x=time,size=nSims,replace=T,prob=q_xpt)*365]

if(sFracType == "CF") {
  
} else if(sFracType == "UDD") {
  events <- events - event_probs[, sample(x=1:365,size=nSims,replace=T)]
}
       

# Sample Lapses               
events2 <- rgeom(nSims,prob=modal_rate2)+1



  


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
  mutate(  Death_Date=Issue_Date %m+% days(events) %m+% days(1),
           Lapse_Date=Issue_Date %m+% months(events2+1) %m+% days(-1)) -> 
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
    .exact_terminations = bExactTerminations
  ) ->
  exposures2

saveRDS(exposures, "exposures.haz.experiment.rds")
exposures <- readRDS("exposures.haz.experiment.rds")

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

census %>%
  inner_join(y=exposures,by=join_by(PolID)) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0),
         Expected_Deaths=-(exposure-2/365)*log(0.5),
         Expected_Deaths2=0.5*(exposure),
         Expected_Deaths3=0.5*(exposure + Death_Count/2),
         Expected_Deaths4=exposure*0.5/(1-exposure*0.5)) %>%
  group_by(Exp_Yr= year(exp_period_end),
           pol_duration) %>%
  summarize(Death_Count=sum(Death_Count),
            Exposure=sum(exposure),
            Expected_Deaths=sum(Expected_Deaths),
            Expected_Deaths2=sum(Expected_Deaths2),
            Expected_Deaths3=sum(Expected_Deaths3),
            Expected_Deaths4=sum(Expected_Deaths4),
            Hazard=Death_Count/sum(exposure-2/365),
            Rate=1-exp(-Hazard),
            Rate2=Death_Count/(Exposure+Death_Count/2),
            A_E_1=Death_Count/Expected_Deaths,
            A_E_2=Death_Count/Expected_Deaths2,
            A_E_3=Death_Count/Expected_Deaths3
  ) %>%
  data.table()

data.table(
  PolID=1:100000,
  Issue_Date=ymd(20200701),
  Death_Date=ymd(20200701) %m+% years(ceiling((events+1)/365)),
  Prem_Mode_Months=12
) -> census_death_padding


census_death_padding %>%
  filter(Issue_Date <= study_end) %>%
  select(PolID,Issue_Date,Death_Date,Prem_Mode_Months) %>%
  expand_exposures(
    .exp_period_start = study_start,
    .exp_period_end = study_end,
    .cal_yr_breaks = cal_yr_breaks,
    .pol_period_granularity = pol_period_granularity,
    .cal_period_granularity = cal_period_granularity,
    .issue_date = Issue_Date,
    .term_date = Death_Date,
    .ID=PolID,
    .prem_mode_months = Prem_Mode_Months
  ) ->
  exposures_death_padding

census_death_padding %>%
  mutate(Death_Date = ymd(20200701) %m+% days(events) %m+% days(1)) %>%
  inner_join(y=exposures,by=join_by(PolID)) %>%
  mutate(Death_Count=ifelse(Death_Date %between% list(exp_period_start,exp_period_end),1,0),
         Expected_Deaths=-exposure*log(0.5),
         Expected_Deaths2=0.5*exposure) %>%
  group_by(Exp_Yr= year(exp_period_end),
           pol_duration) %>%
  summarize(Death_Count=sum(Death_Count),
            Exposure=sum(exposure),
            Expected_Deaths=sum(Expected_Deaths),
            Expected_Deaths2=sum(Expected_Deaths2),
            Hazard=Death_Count/Exposure,
            Rate=1-exp(-Hazard),
            Rate2=Death_Count/(Exposure+Death_Count/2)
  )
