setwd("~/GitHub/fastdid")

library(devtools)
library(tinytest)
library(roxygen2)

load_all()

roxygenise()
build()
run_test_dir()

#before release
build()
check()