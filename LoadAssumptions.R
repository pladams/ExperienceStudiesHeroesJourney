library(data.table)
library(openxlsx)
library(tidyverse)
library(parallel)

rbind(
  read.xlsx(
    xlsxFile = "Assumptions/t3265.xlsx",
    rows=24:102
  ) %>%
    data.table() %>%
    setnames(1,"Issue_Age") %>%
    pivot_longer(
      cols=!Issue_Age,
      names_to="Duration",
      values_to = "qx"
    ) %>%
    mutate(
      Sex="M",
      .before=0
    ) %>%
    filter(!is.na(qx)) %>%
    mutate(
      Issue_Age = as.integer(Issue_Age),
      Duration=as.integer(Duration)
    ) %>%
    data.table(),
  read.xlsx(
    xlsxFile = "Assumptions/t3266.xlsx",
    rows=24:102
  ) %>%
    data.table() %>%
    setnames(1,"Issue_Age") %>%
    pivot_longer(
      cols=!Issue_Age,
      names_to="Duration",
      values_to = "qx"
    ) %>%
    mutate(
      Sex="F",
      .before=0
    ) %>%
    filter(!is.na(qx)) %>%
    mutate(
      Issue_Age = as.integer(Issue_Age),
      Duration=as.integer(Duration)
    ) %>%
    data.table()
) -> vbt2015.select

rbind(
  read.xlsx(
    xlsxFile = "Assumptions/t3265.xlsx",
    rows=116:219
  ) %>%
    setnames(1:2,c("Attained_Age","qx")) %>%
    mutate(Sex="M",.before=0) %>%
    mutate(Attained_Age=as.integer(Attained_Age)),
  read.xlsx(
    xlsxFile = "Assumptions/t3266.xlsx",
    rows=116:219
  ) %>%
    setnames(1:2,c("Attained_Age","qx")) %>%
    mutate(Sex="F",.before=0) %>%
    mutate(Attained_Age=as.integer(Attained_Age))
) -> vbt2015.ult

expand_grid(
 Sex=c("M","F"),
 Issue_Age=18:95,
 Duration=1:120
) %>%
  mutate(Attained_Age = as.integer(Issue_Age + Duration - 1)) %>%
  filter(Attained_Age <= 120) %>%
  left_join(
    y=vbt2015.select
  ) %>%
  rename(qx_sel=qx) %>%
  left_join(
    y=vbt2015.ult
  ) %>%
  rename(qx_ult=qx) %>%
  mutate(qx_sel=as.numeric(qx_sel),
         qx_ult=as.numeric(qx_ult)) %>%
  mutate(qx_su=ifelse(is.na(qx_sel),qx_ult,qx_sel)) %>%
  select(-qx_sel) %>%
  data.table() ->
  vbt2015

vbt2015[Attained_Age==120,
        `:=`(qx_ult=1,qx_su=1)]

rm(vbt2015.select,vbt2015.ult)  

vbt2015[,`:=`(
  tpx_ult=cumprod(1-qx_ult),
  tpx_su=cumprod(1-qx_su)
),
by=.(Sex,Issue_Age)]

vbt2015[
  ,`:=`(
    q_xpt_ult=qx_ult*shift(tpx_ult,fill=1),
    q_xpt_su=qx_su*shift(tpx_su,fill=1)
  ),
  by=.(Sex,Issue_Age)
]

# True mortality
# M 90%
# F 95%
# Duration: 100% for duration 1, then increasing 0.5% per duration for all durations
# Calendar year drift: starting in 2046, deterioration of 1% per year

# Pricing assumption is same, but no duration or calendar adjustment

vbt2015 %>%
  inner_join(
    data.table(
      Duration=1:103,
      Dur_Factor=exp(pmin(20,seq(0,102))*log(1.005))
      )
  ) %>%
  mutate(qx_true=qx_su*Dur_Factor*ifelse(Sex=="M",.9,.95)) %>%
  #select(-qx_ult,-qx_su,-Dur_Factor) %>%
  data.table() ->
  mort.true


sample_death_duration <- function(gender, age, issue_year, mort_factor=1) {
  df <- mort.true[Sex==gender & Issue_Age == age,.(Duration,qx_true)]
  df[,Calendar_Year:=issue_year + Duration - 1]
  df[,qx_true:=pmin(1,mort_factor*qx_true*ifelse(Calendar_Year <= 2045,1,
                              1.01^(pmin(20,Calendar_Year-2045)))
                   )]
  df[age + Duration - 1 ==120,qx_true:=1]
  
  df[,tpx_true:=cumprod(1-qx_true)]
  df[,q_xpt_true:=shift(tpx_true,fill=1)*qx_true]
  
  df[,sample(Duration,1,F,q_xpt_true)]
}

# True lapse
# Very simple:
# High face: Duration 1 is 20%, 2 is 15%, 3-5 is 5%, 6+ is 3%
# Low face: Duration 1 is 10%, 2 is 7%, 3-5 is 5%, 6+ is 3%

rbind(
  data.table(
    Face_Group="High_Face",
    Duration=1:103,
    Lapse_Rate=c(.2,.15,.05,.05,.05,rep(.03,97),1)
  ),
  data.table(
    Face_Group="Low_Face",
    Duration=1:103,
    Lapse_Rate=c(.1,.07,.05,.05,.05,rep(.03,97),1)
  )
) %>%
  arrange(Face_Group,Duration) %>%
  group_by(
    Face_Group
  ) %>%
  mutate(tpx_lapse=cumprod(1-Lapse_Rate)) %>%
  mutate(w_xpt=shift(tpx_lapse,fill=1)*Lapse_Rate) %>%
  data.table() ->
  truelapse



sample_lapse_duration <- function(Face_Group_in) {
  truelapse[Face_Group==Face_Group_in,sample(Duration,1,F,w_xpt)]
}




sample_death_duration_vbt2015_su <- function(gender, age) {
  vbt2015[Sex==gender & Issue_Age==age,sample(Duration,1,F,q_xpt_su)]
}