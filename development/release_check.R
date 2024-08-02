setwd("~/Documents/GitHub/fastdid")

library(devtools)
library(tinytest)
library(roxygen2)

roxygenise()

load_all()

#run test
run_test_dir()

#before release
check()

source("development/build_source.R")
build_source("0.9.4")