rm(list = ls())
gc()

library(devtools)
library(roxygen2)
library(tinytest)
library(usethis)

setwd("~/GitHub/EventStudyCode")

roxygenise()
load_all()
test_all()


# test run