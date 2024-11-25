

library(devtools)
library(tinytest)
library(roxygen2)

roxygenise()

load_all()

#run test
run_test_dir()

#before release
build(path = "development")
check()

source("development/build_source.R")
build_source("1.0.4")