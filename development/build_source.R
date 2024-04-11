rm(list = ls())
gc()

setwd("~/GitHub/fastdid")

source_files <- list.files("R", include.dirs = FALSE, full.names = TRUE)

ver <- "0.9.2"
vername <- "unbalanced"

sink("development/fastdid_sourcever.R")
cat(paste0("#", as.character(Sys.Date()), "\n"))
cat(paste0("message('loading fastdid source ver. ver: ", ver, " (", vername ,"), date: " , as.character(Sys.Date()), "')\n"))
cat("require(data.table);require(stringr);require(BMisc);require(collapse);require(dreamerr);require(parglm)\n")
for(file in source_files){
  current_file = readLines(file)
  cat(current_file, sep ="\n")
  cat("\n")
}

sink()