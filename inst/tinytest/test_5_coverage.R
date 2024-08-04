
#simple coverage
simdt <- sim_did(1e+03, 10, hetero = "all", balanced = TRUE, second_outcome = TRUE)
dt <- simdt$dt
result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time")
cov <- merge(simdt$att, result, by.x = c("G", "time"), by.y = c("cohort", "time"))
cov[, cover := as.numeric(attgt < att_ciub & attgt > att_cilb)]
mean(cov[, cover])

#dynamic coverage
simdt <- sim_did(1e+03, 10, hetero = "dynamic", balanced = TRUE, second_outcome = TRUE)
dt <- simdt$dt
att <- simdt$att
att[, event_time := time - G]
att <- att[, .(attgt = mean(attgt)), by = "event_time"]
result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "dynamic")
cov <- merge(att, result, by = "event_time")
cov[, cover := as.numeric(attgt < att_ciub & attgt > att_cilb)]
mean(cov[, cover])