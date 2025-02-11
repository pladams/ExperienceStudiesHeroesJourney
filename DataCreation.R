library(data.table)
library(plyr)
library(tidyverse)
library(arrow)
library(parallel)
library(rvinecopulib)

RNGkind("L'Ecuyer-CMRG")

source("LoadAssumptions.R")

ilec_path <- "/workspace/Projects/ILEC/VBT/Data/ilecdata_20240119"
policy_census_size <- 2000000
Issue_Year_Range <- 2041:2054

ilec_dataset <- arrow::open_dataset(
  sources=ilec_path,
  format="parquet"
)

ilec_dataset %>%
  filter(Insurance_Plan == "Term" &
           Duration == 1 &
           Issue_Age >= 18 &
           Issue_Age <= 70 &
           Number_of_Pfd_Classes == "4" &
           Observation_Year >= 2011 &
           Observation_Year <= 2017) %>%
  group_by(Sex,Issue_Age,Face_Amount_Band) %>%
  summarize(Policies_Exposed=sum(Policies_Exposed)) %>%
  collect() %>%
  data.table() ->
  src_distribution

src_distribution %>%
  mutate(Face_Group=fct_collapse(Face_Amount_Band,
                                 Low_Face=c("01: 0 - 9,999","02: 10,000 - 24,999","03: 25,000 - 49,999","04: 50,000 - 99,999",
                                            "05: 100,000 - 249,999","06: 250,000 - 499,999"),
                                 other_level="High_Face")
  ) %>%
  group_by(Sex,Issue_Age,Face_Group) %>%
  summarize(Policies_Exposed=sum(Policies_Exposed)) %>%
  arrange(Sex,Issue_Age,Face_Group) %>%
  data.table() ->
  src_distribution

src_distribution %>%
  cross_join(
    data.table( Issue_Year=Issue_Year_Range,
                Proportion=seq(.8,1.3,length=length(Issue_Year_Range)))
  ) %>%
  mutate(Policies_Exposed=Policies_Exposed*Proportion) %>%
  select(-Proportion) ->
  src_distribution

set.seed(0xBEEF)
policy_pop <- src_distribution[
  sample.int(n=nrow(src_distribution),
             size=policy_census_size,
             replace=T,
             prob=src_distribution$Policies_Exposed),
  .(Sex,Issue_Age,Face_Group,Issue_Year)
]

# Issue Date
data.table(
  Month=1:12,
  MaxDays=c(31,28,31,30,31,30,31,31,30,31,30,31),
  Proportion=rep(1/12,12)
) %>%
  cross_join(
    y=data.table(Day=1:31)
  ) %>%
  filter(
    Day<=MaxDays
  ) %>%
  select(-MaxDays) ->
  issue_days

issue_days[Month==1,Proportion:=Proportion*.5]
issue_days[Month==2,Proportion:=Proportion/.7]
issue_days[Month==11,Proportion:=Proportion/.7]
issue_days[Month==12,Proportion:=Proportion/.5]

policy_pop <- cbind(policy_pop,
                          issue_days[
                            sample.int(n=nrow(issue_days), size=policy_census_size,replace=T,prob=issue_days$Proportion),
                            .(Month,Day)
                          ]
)

policy_pop %>%
  mutate(Issue_Date=ymd(paste0(Issue_Year," ",Month," ",Day))) %>%
  select(-Month,-Day) ->
  policy_pop

# Face Amounts
# Low face is 100-450K spaced at 50K intervals
# High Face is 500K - 2M spaced at 100K intervals
policy_pop[Face_Group=="Low_Face",
           Face_Amount:=1000*sample(x=seq(100,450,50),
                               size=nrow(.SD),
                               replace=T)
]

policy_pop[Face_Group=="High_Face",
           Face_Amount:=1000*sample(x=seq(500,2000,100),
                                    size=nrow(.SD),
                                    replace=T)
]


# Create synthetic preferred distributions
# BMI, oil pressure, and oil viscosity are modeled with a copula, which is necessarily
# a cvine (3 variables)
pref.vine <- vinecop_dist(
  pair_copulas = list(
    list(
      bicop_dist("bb7",180,c(1.5,3)),
      bicop_dist("bb7",180,c(1.5,3))
    ),
    list(bicop_dist("gaussian",0,.4))
  ),
  structure = cvine_structure(1:3)
)

rvinecop(n=policy_census_size,
         pref.vine,
         cores=8) ->
  pref.samples



# Variable 3 is BMI, 2 is oil pressure, 1 is oil viscosity


  
pref.samples %>%
  as.data.table() %>%
  rename(
    BMI_u=V3,
    OilViscosity_u=V1,
    OilPressure_u=V2
  ) -> pref.samples

# pref.samples are sampled in the unit cube.

policy_pop %>%
  cbind(pref.samples) ->
  policy_pop





# High Face
# BMI is gamma with mean 25 and variance 50
# plot(seq(10,50,.1),dgamma(seq(10,50,.1),shape=25*25/30,scale=30/25))
# Oil Pressure is gamma with mean 30 and variance 20
# plot(seq(10,50,.1),dgamma(seq(10,50,.1),shape=30*30/20,scale=20/30))
# Oil viscosity is gamma with mean 20 and variance 20
# plot(seq(10,50,.1),dgamma(seq(10,50,.1),shape=20*20/20,scale=20/20))

# Low Face is 90% of these, Males are 120% of these
# Starting in year 2046, apply 1% per issue year increase to all means
low_face_factor_true <- 0.9
male_factor_true <- 1.2

convert_marginal_to_gamma <- function(x, mean=1, variance=1, rounding=2) {
  shape <- mean*mean/variance
  scale <- variance/mean
  
  round(qgamma(x,shape=shape,scale=scale),rounding)
}

policy_pop[,PrefMeanFactor:=ifelse(Sex=="Male",male_factor_true,1)*ifelse(Face_Group=="Low_Face",low_face_factor_true,1)*(1+pmax(0,Issue_Year-2045)/100)]

policy_pop[,
           `:=`(
             BMI=mapply(FUN=convert_marginal_to_gamma,
                        BMI_u,
                        25*PrefMeanFactor,
                        variance=50),
             OilPressure=mapply(FUN=convert_marginal_to_gamma,
                                OilPressure_u,
                                30*PrefMeanFactor,
                                variance=20),
             OilViscosity=mapply(FUN=convert_marginal_to_gamma,
                                OilViscosity_u,
                                20*PrefMeanFactor,
                                variance=20)
           )]


policy_pop[,`:=`(OilViscosity_u=NULL,OilPressure_u=NULL,BMI_u=NULL,PrefMeanFactor=NULL)]

silu <- function(x) {
  x/(1+exp(-x))
}

softplus <- function(x) {
  log(1+exp(x))
}

# Mortality factors for preferreds
# BMI follows a j-curve where every 5 points of BMI increases mortality by 25%
plot(seq(15,35,.1),
     exp(log(1.25)*(seq(15,35,.1)-25)/5)
)

# Oil viscosity and oil pressure follow u-curves modeled as overlaid j-curves
# Oil pressure
x_range <- seq(10,30,.1)


plot(x_range,
     (exp(log(1.25)*softplus( .4*(x_range-20)/1 ))+exp(log(1.05)*softplus( -1*(x_range-20)/1 )))/(1.05+1.25)
)
     

#Oil Viscosity - low viscosity is not as much an issue as low pressure
x_range <- seq(20,40,.1)
plot(x_range,
     (exp(log(1.25)*softplus( .6*(x_range-30)/1 ))+exp(log(1.05)*softplus( -2*(x_range-30)/1 )))/(1.05+1.25)
)


policy_pop[,
           `:=`(PrefFactor_BMI=exp(log(1.25)*(BMI-25)/5),
                PrefFactor_Pressure=(exp(log(1.25)*softplus( .25*(OilPressure-20)/1 ))+exp(log(1.05)*softplus( -1*(OilPressure-20)/1 )))/(1.05+1.25),
                PrefFactor_Viscosity=(exp(log(1.25)*softplus( 0.6*(OilViscosity-30)/1 ))+exp(log(1.05)*softplus( -2*(OilViscosity-30)/1 )))/(1.05+1.25)
                )]

# The softplus unfortunately isn't automatically normalized
policy_pop[,
           `:=`(PrefFactor_Viscosity=PrefFactor_Viscosity/mean(PrefFactor_Viscosity),
                PrefFactor_Pressure=PrefFactor_Pressure/mean(PrefFactor_Pressure))]

policy_pop[,PrefFactor:=PrefFactor_BMI*PrefFactor_Pressure*PrefFactor_Viscosity]

# Preferred criteria finding - traditional
policy_pop[,quantile(BMI,.95)]
policy_pop[,quantile(OilPressure,.95)]
policy_pop[,quantile(OilViscosity,.95)]

# Definition of standard
policy_pop[BMI <= 41 & OilPressure <= 41 & OilViscosity <= 32,.(sum(Face_Amount*PrefFactor)/sum(Face_Amount),sum(Face_Amount),.N)]

policy_pop[,
           UW_Decision:="SUB"]
policy_pop[BMI <= 41 & OilPressure <= 41 & OilViscosity <= 32,
           UW_Decision:="STD"]
policy_pop[UW_Decision == "SUB" & PrefFactor > 3, 
           UW_Decision := "DEC"]

policy_pop[,PrefClass:=3]
policy_pop[BMI <= 30 & OilPressure <= 35 & OilViscosity <= 25 & UW_Decision == "STD",
           PrefClass:=1]
policy_pop[PrefClass==1,
           .(sum(Face_Amount*PrefFactor)/sum(Face_Amount),sum(Face_Amount),.N)]

policy_pop[BMI <= 33 & OilPressure <= 38 & OilViscosity <= 28 & PrefClass > 1 & UW_Decision == "STD",
           PrefClass:=2]
policy_pop[PrefClass == 2 & BMI <= 33 & OilPressure <= 38 & OilViscosity <= 28,.(sum(Face_Amount*PrefFactor)/sum(Face_Amount),sum(Face_Amount),.N)]


policy_pop[PrefClass == 3 & UW_Decision == "STD",.(sum(Face_Amount*PrefFactor)/sum(Face_Amount),sum(Face_Amount),.N)]


policy_pop[UW_Decision == "STD" & year(Issue_Date) < 2046,.(sum(Face_Amount*PrefFactor)/sum(Face_Amount),sum(Face_Amount),.N),by=.(PrefClass)][order(PrefClass)]



policy_pop[,
           Death_Duration:=mcmapply(FUN=sample_death_duration,
                  Sex,
                  Issue_Age,
                  Issue_Date,
                  PrefFactor,
                  mc.cores = 16)]


policy_pop[,
           Lapse_Duration:=mcmapply(FUN=sample_lapse_duration,
                                    Face_Group,
                                    mc.cores = 16)]

policy_pop[,prem_mode:= sample(c("A","Q","M"),size=nrow(.SD),prob=c(.04,.01,.95),replace=T)]

# Cutoff 12-31-2053
death.skew <- 1:365
death.skew <- (1+death.skew*.1/364-.05-.1/364)/365

policy_pop[,death_skew_days:= sample.int(365,size=nrow(.SD),replace=T,prob = death.skew)]
policy_pop[,lapse_skew_months := sample.int(12,size=nrow(.SD),replace=T)]

policy_pop[,
           Death_Date:=Issue_Date %m+% years(Death_Duration-1) %m+% days(death_skew_days-1)]

policy_pop[prem_mode=="M",
           Lapse_Date:=Issue_Date %m+% years(Lapse_Duration-1) %m+% months(lapse_skew_months)]
policy_pop[prem_mode=="A",
           Lapse_Date:=Issue_Date %m+% years(Lapse_Duration-1)]
policy_pop[prem_mode=="Q",
           Lapse_Date:=Issue_Date %m+% years(Lapse_Duration-1) %m+% months(3*(lapse_skew_months %% 4 + 1 ) )]

policy_pop[,Status:="Active"]
policy_pop[Death_Date <= ymd(20531231) & Death_Date < Lapse_Date, `:=`(Status="Death",Term_Date=Death_Date)]
policy_pop[Lapse_Date <= ymd(20531231) & Lapse_Date < Death_Date, `:=`(Status="Lapsed",Term_Date=Lapse_Date)]

policy_pop[UW_Decision=="DEC",
           `:=`(
             Status="Not Issued",
             Term_Date=NA
           )]

policy_pop[,PolID:=1:nrow(.SD)]






# Showing that hazard method recovers qx
policy_pop[,`:=`(Sex_HazTest="M",Issue_Age_HazTest=65)]

Death_Duration_VBT2015_samp <- vbt2015[Sex=="M" & Issue_Age==65,sample(x=Duration,size=nrow(policy_pop),replace=T,prob=q_xpt_su)]

policy_pop[,
           Death_Duration_VBT2015:=Death_Duration_VBT2015_samp]
policy_pop[,
           Death_Date_VBT2015SU:=Issue_Date %m+% years(Death_Duration_VBT2015-1) %m+% days(death_skew_days-1)]

policy_pop[,Status_ValHazardMethod:="Active"]
policy_pop[Death_Date_VBT2015SU <= ymd(20531231) & Death_Date_VBT2015SU < Lapse_Date, `:=`(Status_ValHazardMethod="Death",Term_Date_ValHazMethod=Death_Date_VBT2015SU)]
policy_pop[Lapse_Date <= ymd(20531231) & Lapse_Date < Death_Date_VBT2015SU, `:=`(Status_ValHazardMethod="Lapsed",Term_Date_ValHazMethod=Lapse_Date)]


# policy_pop[,
#            `:=`(
#              Death_Duration=NULL,
#              Lapse_Duration=NULL,
#              Death_Date=NULL,
#              Lapse_Date=NULL,
#              Death_Date_VBT2015SU=NULL
#            )]


saveRDS(policy_pop, "policy_pop.rds")
write_csv(
  x=policy_pop,
  file="policy_pop.csv",
  col_names = F,
  na=""
)
