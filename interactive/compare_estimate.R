
rm(list = ls())
gc()

library(devtools)
library(profvis)

setwd("~/GitHub/EventStudyCode")

load_all()

library(did)

# setup --------------------------------

sample_size <- 1000
time_period <- 10
dt <- sim_did(sample_size, time_period, seed = 1, hetero = "dynamic")[["dt"]]


#cohort event time est ----------------------

event_panel <- suppressWarnings(create_event_data(dt, timevar = "time", unitvar =  "unit", cohortvar = "G",
                                                  control_group = "both", copy = TRUE, verbose = FALSE))
event_est_ce <- get_event_result(event_panel, variable = "y", trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")
event_est_ce[, method := "ec"]

cohort_obs <- dt[!is.infinite(G), .(count = uniqueN(unit)) , by = "G"]
event_est_ce <- event_est_ce |> merge(cohort_obs, by.x = "cohort", by.y = "G")

event_est_sum <- event_est_ce[, .(Estimate = sum(Estimate*count)/sum(count),
                                  `Std. Error` = sqrt(sum(`Std. Error`^2*count*2))/sum(count)) , by = "event_time"]
event_est_sum[, method := "ec_aggregated"]

#dynamic est ----------------------
event_panel <- suppressWarnings(create_event_data(dt, timevar = "time", unitvar =  "unit", cohortvar = "G",
                                                  control_group = "both", copy = TRUE, verbose = FALSE))
event_est_dynamic <- get_event_result(event_panel, variable = "y", trends = FALSE, mem.clean = FALSE, result_type = "dynamic")
event_est_dynamic[, method := "ec_dynamic"]

#did estimates ------------------
did_result <- att_gt(yname = "y",
                  gname = "G",
                  idname = "unit",
                  tname = "time",
                  data = dt,
                  base_period = "universal",
                  control_group = "notyettreated")
did_est <- data.table(Estimate = did_result$att, `Std. Error` = did_result$se, cohort = unlist(did_result$group), time = unlist(did_result$t))
did_est[, event_time := time - cohort]
did_est[, method := "did"]
did_est[, time := NULL]

#dynamic 
did_dynamic_result <- aggte(did_result, "dynamic", clustervars = "unit")
did_dynamic_est <- data.table(Estimate = did_dynamic_result$att.egt , `Std. Error` = did_dynamic_result$se.egt, event_time = did_dynamic_result$egt)
did_dynamic_est[, method := "did"]

#compare ---------------------------------------

est_compare <- rbind(did_est, event_est_ce[,.(cohort, event_time, Estimate, `Std. Error`, method)])
est_compare |> ggplot(aes(x = method, y = Estimate)) + geom_point() + geom_errorbar(aes(ymax = Estimate + 1.96*`Std. Error`,
                                                                                            ymin = Estimate - 1.96*`Std. Error`)) + 
  facet_wrap(~event_time + cohort, scales = "free")


dynamic_compare <-  rbind(did_dynamic_est, event_est_dynamic[,.(event_time, Estimate, `Std. Error`, method)],
                          event_est_sum[,.(event_time, Estimate, `Std. Error`, method)])

dynamic_compare |> ggplot(aes(x = method, y = Estimate)) + geom_point() + geom_errorbar(aes(ymax = Estimate + 1.96*`Std. Error`,
                                                                                                ymin = Estimate - 1.96*`Std. Error`)) + 
  facet_wrap(~event_time, scales = "free")

# old event code


source("old/source_raw/eventcode_IV.R")

event_panel <- suppressWarnings(create_event_data_old(dt, timevar = "time", unitvar =  "unit", 
                                                      cohortvar = "G",
                                                      balanced_panel = FALSE,
                                                      never_treat_action = "both"))
event_panel <- construct_event_variables(event_panel)
event_es <- get_result_dynamic(event_panel, variable = "y", trends = FALSE)

