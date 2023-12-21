profvis(run_fastdid(dt))

profvis(run_did(dt))

profvis({
  result <- fastdid(dt, 
                    timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                    covariatesvar = "x")
})

#covariates multiple outcome
profvis({
  result <- fastdid(dt, 
                    timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = c("y", "y2"),  result_type = "group_time",
                    covariatesvar = "x")
})

bm_time <- microbenchmark(result <- fastdid(dt, 
                                            timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                                            covariatesvar = "x"),
                          times = 1)
message(bm_time$time/10^9)