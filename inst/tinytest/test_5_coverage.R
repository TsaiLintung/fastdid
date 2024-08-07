
# #gt coverage, unconfounded
# simdt <- sim_did(1e+04, 10, hetero = "all", balanced = TRUE, second_outcome = TRUE)
# dt <- simdt$dt
# result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time")
# cov <- merge(simdt$att, result, by.x = c("G", "time"), by.y = c("cohort", "time"))
# cov[, cover := as.numeric(attgt < att_ciub & attgt > att_cilb)]
# mean(cov[, cover])
# 
# #gt coverage, confounded
# simdt_c <- sim_did(1e+04, 10, hetero = "all", balanced = TRUE, second_outcome = TRUE, second_cohort = TRUE)
# dt_c <- simdt_c$dt
# result <- fastdid(dt_c, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time")
# cov <- merge(simdt_c$att[event == 1], result, by.x = c("G", "time"), by.y = c("cohort", "time"))
# cov[, cover := as.numeric(attgt < att_ciub & attgt > att_cilb)]
# mean(cov[, cover])

#gt coverage, confounded but fixed
covers <- c()
pb <- progress::progress_bar$new(total = 100)
for(i in 1:100){
  simdt_c <- sim_did(1e+03, 10, hetero = "all", balanced = TRUE, second_outcome = TRUE, second_cohort = TRUE)
  dt_c <- simdt_c$dt
  result <- fastdid(dt_c, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                    exper = list(cohortvar2 = "G2", event_specific = TRUE), alpha = 0.05, boot = TRUE, cband = FALSE)
  cov <- merge(simdt_c$att[event == 1], result, by.x = c("G", "time"), by.y = c("cohort", "time"))
  cov[, cover := as.numeric(attgt < att_ciub & attgt > att_cilb)]
  covers <- c(covers, cov[, cover])
  message(i)
}
mean(covers)


#gt coverage, normal
covers <- c()
pb <- progress::progress_bar$new(total = 100)
for(i in 1:100){
  simdt_c <- sim_did(1e+03, 10, hetero = "all", balanced = TRUE, second_outcome = TRUE, second_cohort = FALSE)
  dt_c <- simdt_c$dt
  result <- fastdid(dt_c, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"
                    , alpha = 0.05, boot = TRUE, cband = FALSE)
  cov <- merge(simdt_c$att, result, by.x = c("G", "time"), by.y = c("cohort", "time"))
  cov[, cover := as.numeric(attgt < att_ciub & attgt > att_cilb)]
  covers <- c(covers, cov[, cover])
  pb$tick()
}
mean(covers)