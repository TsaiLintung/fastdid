
rm(list = ls())
gc()

library(profvis)
library(microbenchmark)
library(tinytest)

setwd("~/GitHub/EventStudyCode")

# load event code ---------------------------------------------------------------------

source("sim_did.R")

source("source/setup.R")
source("source/create_event_data.R")
source("source/get_event_result.R")

source("source/test_funcs.R")

# setup --------------------------------------------------------------------------



# test  ---------------------------------------------------------------------

p <- list()

p$time_period <- 10
p$sample_size <- 100
p$min_time <- -Inf
p$max_time <- Inf
p$y_name <- c("y")
p$t_name <- "time"
p$unit_name <- "unit"
p$cohort_name <- "G"
p$stratify_name <- "s"
p$balance_name <- "x"

test_create_event_data()
#test_get_result

# test estimations

p$sample_size <- 1000

test_get_result_dynamic(p)
test_get_result_cohort_event_time(p)
test_get_result_means()