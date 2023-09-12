
rm(list = ls())
gc()

library(profvis)
library(microbenchmark)
library(testthat)

setwd("~/GitHub/EventStudyCode")

# load event code ---------------------------------------------------------------------

source("sim_did.R")

source("source/setup.R")
source("source/preprocess.R")
source("source/estimation.R")
source("source/report.R")
source("source/test.R")

# simulation ---------------------------------------------------------------------


time_period <- 10
sample_size <- 1000

simdt <- sim_did(sample_size, time_period, cov = "int", hetero = "dynamic", balanced = FALSE, second_outcome = FALSE, seed = 1)
dt <- simdt$dt
att <- simdt$att

  
event_panel <- copy(dt) #copying so that the original does not change

min_time <- -Inf
max_time <- Inf
y_name <- c("y")
t_name <- "time"
unit_name <- "unit"
cohort_name <- "G"
stratify_name <- "s"
balance_name <- "x"

event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                                 cohortvar = cohort_name,
                                                 covariate_base_balance = balance_name,
                                                 covariate_base_stratify = stratify_name,
                                                 balanced_panel = TRUE,
                                                 never_treat_action = "both")

dynamic_est <- get_result_dynamic(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)

ce <- get_result_cohort_event_time(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)


pooled_est <- get_result_pooled(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)
means_est <- get_result_means(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)
#event_est <- get_result_covariates(event_panel, covariate = stratify_name, variable = y_name, trends = FALSE, mem.clean = FALSE)

#test for dynamic estimate

test_dynamic_est <- function(dynamic_est, att){
  att[, event_time := time-G]
  att_dynamic <- att[!event_time %in% c(-1,time_period-1), .(attgt = mean(attgt)), by = "event_time"]
  event_validate <- merge(dynamic_est, att_dynamic, by = "event_time")
  event_validate[, ci_ub := Estimate+`Std. Error`*1.96]
  event_validate[, ci_lb := Estimate-`Std. Error`*1.96]
  event_validate[, s :=  (function(x){str_sub(str_trim(x), str_length(str_trim(x)), str_length(str_trim(x)))})(variable)]
  event_validate[, attgt := attgt*as.numeric(s)]
  event_validate[, par_in_ci := (attgt <= ci_ub & attgt >= ci_lb)]
  expect_gt(event_validate[, mean(par_in_ci)], 0.95)
}

test_that("dynamic estimate",
          {test_dynamic_est(dynamic_est, att)})

