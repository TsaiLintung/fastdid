#by LinTung Tsai to be used with Max's event code
#some wrappers
library(fixest)

estimate_event_dynamics <- function(panel, start, end, outcomes, control = c(), stratify = c(), use_never_treat = FALSE,
                               timevar = "time", unitvar = "id", cohortvar = "event_time"){
  
  event_panel <- copy(panel) #copying so that the original does not change
  
  setnames(event_panel, cohortvar, "event_time")
  onset <- event_panel[, min(event_time) - 1]
  event_panel[, onset_time := onset]
  event_panel <- event_panel %>% create_event_data(timevar = timevar, unitvar = unitvar, 
                                             cohortvar = "event_time",
                                             onset_agevar = "onset_time",
                                             covariate_base_balance = control,
                                             covariate_base_stratify = stratify,
                                             never_treat_action = ifelse(use_never_treat, "both", "exclude"),
                                             balanced_panel = TRUE)
  
  event_panel <- event_ATTs_head(event_panel, outcomes)
  
  message("preprocessing done.")
  
  dt_dynamic <- data.table()
  
  for(outcome in outcomes){
    dt_dynamic_outcome <-data.table()
    
    dt_dynamic_outcome <- suppressMessages(get_result_dynamic(event_panel,start,-2,outcome,dt_dynamic_outcome)) #add new result to dt_dynamic
    dt_dynamic_outcome <- suppressMessages(get_result_dynamic(event_panel,0,end,outcome,dt_dynamic_outcome)) #add new result to dt_dynamic
    
    dt_dynamic_outcome[, outcome_type := outcome]
    dt_dynamic <- rbind(dt_dynamic, dt_dynamic_outcome)
    message(outcome, " done.")
  }
  
  return(dt_dynamic)
}

plot_event_study <-function(dt, graphname, note = ""){
  
  # Author : 冠儒, Hsiao, Don
  # Motified from graph 2 and graph2 0720
  
  #add the base period
  for (outcome in dt[, unique(outcome_type)]) {
    dt <- add_row(dt, variable = "treated_event_time_stratify-1.0", 
                  model = 1, Estimate = 0, "Std. Error" = 0, "t value" = 0, "Pr(>|t|)" = 0, obs = 196890,
                  outcome_type = outcome)
    dt <- add_row(dt, variable = "treated_event_time_stratify-1.1", 
                  model = 1, Estimate = 0, "Std. Error" = 0, "t value" = 0, "Pr(>|t|)" = 0, obs = 196890,
                  outcome_type = outcome)
  }
  
  significance_level <- 0.05
  dt[, c("first", "stratify_dummy") := tstrsplit(variable, ".", fixed = TRUE)]
  dt[, coefficient := str_extract(first, '([a-z\\_]+)')]
  dt[, event_time := str_extract(first, '([-0-9]+)')]
  dt[, event_time := as.numeric(event_time)]
  dt[, conf_upb := Estimate + qnorm(1-significance_level/2) * `Std. Error`]
  dt[, conf_lwb := Estimate - qnorm(1-significance_level/2) * `Std. Error`]
  
  figure <- dt[str_sub(variable, 1, 7) == "treated"] %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'red') +
    geom_line(aes(x = event_time, y = Estimate, color = "black")) + 
    geom_point(aes(x = event_time, y = Estimate, color = "black")) +
    geom_errorbar(aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
                  width = 0.1, linetype = "dashed") +
    facet_wrap(~ outcome_type, scales = "free") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.background = element_rect(linetype = "dashed", color = "black"),
          legend.box = "horizontal",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0)) +
    #scale_y_continuous(name = ylabel, limits = yscale, breaks = ybks) +
    #scale_x_continuous(name = xlabel, limits = xscale, breaks = xbks) +
    labs(title = graphname, subtitle = note)

  return(figure)
  
}





