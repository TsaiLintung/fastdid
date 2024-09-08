#setup

#takes a lot of time so only at home
if(at_home() & .Platform$OS.type == "unix" & requireNamespace("parallel")){
  
  #gotta go fast 
  setDTthreads(0)
  options(mc.cores = getDTthreads())
  
  
  cov_error_margin <- 0.01
  iter <- 200
  
  #doubledid boot
  covers <- c()
  pb <- progress::progress_bar$new(total = iter)
  for(i in 1:iter){
    simdt_c <- sim_did(1e+04, 5, hetero = "all", balanced = TRUE, second_outcome = TRUE, second_cohort = TRUE)
    dt_c <- simdt_c$dt
    result <- fastdid(dt_c, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", event_specific = TRUE, alpha = 0.05, boot = TRUE, cband = FALSE, parallel = TRUE)
    cov <- merge(simdt_c$att[event == 1], result, by.x = c("G", "time"), by.y = c("cohort", "time"))
    cov[, cover := as.numeric(attgt < att_ciub & attgt > att_cilb)]
    covers <- c(covers, cov[, cover])
    pb$tick()
  }
  mc <- mean(covers)
  # the standard error is sqrt(0.95*0.05/3200), so may wrongfully reject about 1% of time
  expect_true(mc > 0.95-cov_error_margin & mc < 0.95+cov_error_margin, info = "double did, bootstrap coverage")
  
  rm(mc, covers, simdt_c, dt_c, cov)
  
  #double did, analytical 
  pb <- progress::progress_bar$new(total = iter)
  covers <- c()
  for(i in 1:iter){
    simdt_c <- sim_did(1e+04, 5, hetero = "all", balanced = TRUE, second_outcome = TRUE, second_cohort = TRUE)
    dt_c <- simdt_c$dt
    result <- fastdid(dt_c, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", event_specific = TRUE, alpha = 0.05, boot = FALSE, cband = FALSE, parallel = TRUE)
    cov <- merge(simdt_c$att[event == 1], result, by.x = c("G", "time"), by.y = c("cohort", "time"))
    cov[, cover := as.numeric(attgt < att_ciub & attgt > att_cilb)]
    covers <- c(covers, cov[, cover])
    pb$tick()
  }
  mc <- mean(covers)
  # the standard error is sqrt(0.95*0.05/3200), so may wrongfully reject about 1% of time
  expect_true(mc > 0.95-cov_error_margin & mc < 0.95+cov_error_margin, info = "double did, analytical coverage")
  
  #unbalanced panel OR
  pb <- progress::progress_bar$new(total = iter)
  covers <- c()
  for(i in 1:iter){
    simdt_c <- sim_did(1e+03, 10, hetero = "all", balanced = TRUE, second_outcome = TRUE, cov = "cont")
    dt <- simdt_c$dt
    dt2 <- data.table::copy(dt)
    keep <- sample(c(rep(TRUE, 15),FALSE), dt2[,.N], TRUE)
    dt2 <- dt2[keep]
    result <- fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      allow_unbalance_panel = TRUE, covariatesvar = "x", control_type = "reg", parallel = TRUE)
    cov <- merge(simdt_c$att, result, by.x = c("G", "time"), by.y = c("cohort", "time"))
    cov[, cover := as.numeric(attgt < att_ciub & attgt > att_cilb)]
    covers <- c(covers, cov[, cover])
    pb$tick()
  }
  mc <- mean(covers)
  # the standard error is sqrt(0.95*0.05/3200), so may wrongfully reject about 1% of time
  expect_true(mc > 0.95-cov_error_margin & mc < 0.95+cov_error_margin, info = "unbalanced panel, regression")
}
