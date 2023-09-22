
sim_dt <- sim_did(1000, 10, hetero = "dynamic", seed = 1) #generate dataset
att <- sim_dt$att
dt <- sim_dt$dt
event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit",
                                 cohortvar = "G",
                                 control_group = "both", verbose = FALSE)
dynamic_est <- event_panel |> get_event_result(variable = "y", result_type = "dynamic") #estimate
att[, event_time := time-G] #process att 
att_dynamic <- att[!event_time %in% c(-1,10-1), .(attgt = mean(attgt)), by = "event_time"]
merged_att <- merge(dynamic_est, att_dynamic, by = "event_time")
ratio <- get_att_in_ci_ratio(merged_att)
expect_true(ratio>0.9, info = "(almost all) att in dynamic's CI, no balance and stratify")

cohort_time_est <- event_panel |> get_event_result(variable = "y", result_type = "cohort_event_time") #estimate
att[, event_time := time-G]
cohort_time_att <- att[!event_time %in% c(-1,10-1),]
get_cohort <- function(x){
  dot_pos <- str_locate_all(x, "\\.")
  start <- dot_pos[[1]][1, 1]
  end <- dot_pos[[1]][2, 1]
  return(as.numeric(str_sub(x, start + 1, end - 1)))
}
cohort_time_est[,G:= lapply(variable, get_cohort)]
cohort_time_est[, G:= unlist(G)]
merged_att <- merge(cohort_time_est, cohort_time_att, by = c("event_time", "G"))
ratio <- get_att_in_ci_ratio(merged_att)
expect_true(ratio>0.9, info = "(almost all) att in cohort_event_time's CI, no balance and stratify")

expect_silent(event_panel |> get_event_result(variable = "y", result_type = "pooled"),
              info = "pooled result run without error, no balance and stratify")
expect_silent(event_panel |> get_event_result(variable = "y", result_type = "means"),
              info = "mean result run without error, no balance and stratify")
expect_silent(event_panel |> get_event_result(variable = "y",covariate = "x", result_type = "covariate"),
              info = "covariate result run without error, no balance and stratify")


rm(sim_dt, att, dt, event_panel, dynamic_est, att_dynamic, merged_att, ratio, cohort_time_att, cohort_time_est)


sim_dt <- sim_did(1000, 10, cov = "int", stratify = TRUE, hetero = "dynamic", seed = 1) #generate dataset
att <- sim_dt$att
dt <- sim_dt$dt
event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit",
                                 cohortvar = "G",
                                 covariate_base_stratify = "s",
                                 covariate_base_balance = "x",
                                 control_group = "both", verbose = FALSE)
dynamic_est <- event_panel |> get_event_result(variable = "y", result_type = "dynamic") #estimate
att[, event_time := time-G] #process att 
att_dynamic <- att[!event_time %in% c(-1,10-1), .(attgt = mean(attgt)), by = "event_time"]
merged_att <- merge(dynamic_est, att_dynamic, by = "event_time")
ratio <- get_att_in_ci_ratio(merged_att)
expect_equal(ratio, 1, info = "att in dynamic's CI, with balance and stratify")

cohort_time_est <- event_panel |> get_event_result(variable = "y", result_type = "cohort_event_time") #estimate
att[, event_time := time-G]
cohort_time_att <- att[!event_time %in% c(-1,10-1),]
get_cohort <- function(x){
  dot_pos <- str_locate_all(x, "\\.")
  start <- dot_pos[[1]][1, 1]
  end <- dot_pos[[1]][2, 1]
  return(as.numeric(str_sub(x, start + 1, end - 1)))
}
cohort_time_est[,G:= lapply(variable, get_cohort)]
cohort_time_est[, G:= unlist(G)]
merged_att <- merge(cohort_time_est, cohort_time_att, by = c("event_time", "G"))
ratio <- get_att_in_ci_ratio(merged_att)
expect_true(ratio>0.9, info = "(almost all) att in cohort_event_time's CI, with balance and stratify")

expect_silent(event_panel |> get_event_result(variable = "y", result_type = "dynamic", trends = TRUE),
              info = "dynamic result run without error, w balance and stratify, add trends")
expect_silent(event_panel |> get_event_result(variable = "y", result_type = "cohort_event_time", trends = TRUE),
              info = "cohort_event_time result run without error, w balance and stratify, add trends")
expect_silent(event_panel |> get_event_result(variable = "y", result_type = "pooled"),
              info = "pooled result run without error, w balance and stratify")
expect_silent(event_panel |> get_event_result(variable = "y", result_type = "means"),
              info = "mean result run without error, w balance and stratify")
expect_silent(event_panel |> get_event_result(variable = "y",covariate = "x", result_type = "covariate"),
              info = "covariate result run without error, w balance and stratify")

rm(sim_dt, att, dt, event_panel, dynamic_est, att_dynamic, merged_att, ratio, cohort_time_att, cohort_time_est, get_cohort)

