expand_exposures <- function(.data,
                             .exp_period_start,
                             .exp_period_end,
                             .cal_yr_breaks,
                             .pol_period_granularity,
                             .cal_period_granularity,
                             .issue_date,
                             .term_date,
                             .ID,
                             .prem_mode_months) {
  
  .data %>%
    mutate(
      first_date=as.Date(
        pmax({{.issue_date}},.exp_period_start)
      ),
      last_date=as.Date(
        pmin(replace_na({{.term_date}},.exp_period_end),.exp_period_end)
      ),
      full_pol_periods = interval({{.issue_date}},last_date) / months(.pol_period_granularity),
      prestudy_pol_periods = interval({{.issue_date}},first_date) / months(.pol_period_granularity)
    ) ->
    .data
  
  modal_sequence <- pol_period_granularity*as.numeric(seq(to=.data[,ceiling(max(full_pol_periods))]))
  
  .data %>%
    select({{.ID}},
           prestudy_pol_periods,
           full_pol_periods,
           {{.prem_mode_months}},
           {{.issue_date}},
           first_date) %>%
    cross_join(y=data.table(monthaversary=modal_sequence)) %>%
    filter(monthaversary >= ceiling(prestudy_pol_periods*.pol_period_granularity) &
             monthaversary <= floor(full_pol_periods*.pol_period_granularity)) %>%
    filter(monthaversary %% {{.prem_mode_months}} == 0) %>%
    mutate(exp_period_end = {{.issue_date}} %m+% months(monthaversary) %m+% days(-1)) %>%
    select({{.ID}},first_date,monthaversary,exp_period_end,{{.issue_date}}) ->
    stage1a
  
  .data %>%
    filter(!is.na({{.term_date}})) %>%
    filter({{.term_date}} >= study_start & {{.term_date}} <= study_end) %>%
    select({{.ID}},{{.issue_date}},first_date,{{.term_date}}) %>%
    rename(exp_period_end={{.term_date}} ) %>%
    mutate(monthaversary=interval({{.issue_date}},exp_period_end) / months(1)) %>%
    select({{.ID}},first_date,monthaversary,exp_period_end,{{.issue_date}}) ->
    stage1b
  
  .data %>%
    cross_join(data.table(exp_period_end=.cal_yr_breaks)) %>%
    filter(exp_period_end >= first_date &
             exp_period_end <= last_date) %>%
    mutate(monthaversary=interval({{.issue_date}},exp_period_end) / months(1)) %>%
    select({{.ID}},first_date,monthaversary,exp_period_end,{{.issue_date}}) ->
    stage1c
  
  rbind(
    stage1a,
    stage1b,
    stage1c
  ) %>%
    distinct() %>%
    arrange({{.ID}}, exp_period_end) %>%
    group_by({{.ID}}) %>%
    mutate(exp_period_start=lag(exp_period_end) %m+% days(1),
           .before=exp_period_end) %>%
    mutate(exp_period_start=as.Date(ifelse(is.na(exp_period_start),first_date,exp_period_start))) %>%
    mutate(pol_duration = pmax(1,ceiling(interval({{.issue_date}},exp_period_end) / years(1))),
           exposure = (interval(exp_period_start,exp_period_end) / days(1) + 1) / ifelse(year(exp_period_end) %% 4 == 0, 366, 365)
    ) %>%
    select({{.ID}},monthaversary,exp_period_start,exp_period_end,pol_duration,exposure) %>%
    data.table() ->
    .data
  
  .data
}