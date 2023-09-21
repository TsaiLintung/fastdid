library(devtools)
library(roxygen2)
library(tinytest)
library(usethis)




setwd("~/GitHub/EventStudyCode")



roxygenise()
load_all()

run_test_file("inst/tinytest/test_create_event_data.R", verbose=0)

