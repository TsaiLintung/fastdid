
rm(list = ls())
gc()

library(profvis)

setwd("~/GitHub/EventStudyCode")

# load event code ---------------------------------------------------------------------

source("sim_did.R")
source("source/setup.R")


simdt <- sim_did(100000, 10, cov = "int", hetero = "dynamic", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = FALSE)
dt <- simdt$dt

min_time <- -Inf			
max_time <- Inf			
y_name <- "y"			
t_name <- "time"			
unit_name <- "unit"			
cohort_name <- "G"			
balance_covariate <- "x"
balance_name <- "x"

# old setup -------------------------------------------------------------------------

source("source/oldfunc.R")

event_panel <- copy(dt) #copying so that the original does not change			

profvis({
  event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 			
                                                   cohortvar = cohort_name,			
                                                   covariate_base_balance = balance_covariate,			
                                                   balanced_panel = TRUE,			
                                                   never_treat_action = "both")		 
})

profvis({			
  event_est <- get_result_dynamic(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)			
})



profvis({			
  event_est <- get_event_result(event_panel, variable = y_name, result_type = "dynamic")			
})

# testing -------------------------------------------------------------

source("source/create_event_data.R")

event_panel <- copy(dt) #copying so that the original does not change			


profvis({
  event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                                   cohortvar = cohort_name,
                                                   covariate_base_balance = balance_name,
                                                   balanced_panel = TRUE,
                                                   control_group = "both")
})

contruct_event_var <- function(eventdata){
  
  base_time <- -1
  
  eventdata[,unitfe := finteraction(time_pair,treated,stratify)]
  
  base_event_stratify <- paste0(c(base_time,1),collapse=".")
  
  eventdata[,event_time_stratify:= finteraction(event_time_fact,stratify)]
  eventdata[,treated_event_time_stratify := event_time_stratify]
  
  #Omitting base year for all levels of --stratify--:
  eventdata[event_time==base_time, event_time_stratify := base_event_stratify]
  eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = base_event_stratify)]
  
  #Omitting base year for all levels of --stratify--, for treated people
  eventdata[event_time==base_time | treated == 0 ,treated_event_time_stratify := base_event_stratify]
  eventdata[,treated_event_time_stratify:=relevel(treated_event_time_stratify,ref = base_event_stratify)]
  
}


