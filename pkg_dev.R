library(devtools)
library(roxygen2)
library(tinytest)
library(usethis)

require(data.table)
require(stringr)
require(fixest)
require(rlist)
require(collapse)
require(dreamerr)



setwd("~/GitHub/EventStudyCode")

load_all()

roxygenise()

run_test_file("inst/tinytest/test_plot_event_dynamics.R", verbose=0)

