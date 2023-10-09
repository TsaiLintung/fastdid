rm(list = ls())
gc()

library(devtools)
library(roxygen2)
library(tinytest)
library(usethis)

setwd("~/GitHub/EventStudyCode")

if(FALSE){
  roxygenise()
  load_all()
  test_all()
}


if(FALSE){
  use_github_actions()
  use_github_action_check_standard()
}
# test run