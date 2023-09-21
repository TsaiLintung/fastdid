
sim_dt <- sim_did(1000, 10, hetero = "dynamic", seed = 1)
att <- sim_dt$att
dt <- sim_dt$dt

event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                 cohortvar = "G",
                                 control_group = "both", verbose = FALSE)
dynamic_est <- event_panel |> get_event_result(variable = "y", result_type = "dynamic")

#process att 
att[, event_time := time-G]
att_dynamic <- att[!event_time %in% c(-1,10-1), .(attgt = mean(attgt)), by = "event_time"]

#merge both
merged_att <- merge(dynamic_est, att_dynamic, by = "event_time")

ratio <- get_att_in_ci_ratio(merged_att)

expect_equal(ratio, 1, info = "att in dynamic's CI, no covariate and stratify")

rm(sim_dt, att, dt, event_panel, dynamic_est, att_dynamic, merged_att, ratio)