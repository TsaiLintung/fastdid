setwd("~/Documents/GitHub/fastdid")

library(devtools)
library(tinytest)
library(roxygen2)

roxygenise()

load_all()
run_test_dir()

#before release
check()