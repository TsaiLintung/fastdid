rm(list = ls())
gc()

library(stringr)
library(here)

source_files <- list.files(here("R"), include.dirs = FALSE, full.names = TRUE)

ver <- "0.9.2"
vername <- "unbalanced"
dep <- c("data.table", "stringr", "BMisc", "collapse", "dreamerr", "parglm")

sink(here("development/source_head.R"))
cat(paste0("#", as.character(Sys.Date()), "\n"))
cat(paste0("message('loading fastdid source ver. ver: ", ver, " (", vername ,"), date: " , as.character(Sys.Date()), "')\n"))
cat(paste0("require(", dep, ");") |> str_flatten())
sink()

sink(paste0("development/fastdid_", str_replace_all(ver, "\\.", "_"), ".R"))
cat(readLines(here("development/source_head.R")), sep ="\n")
for(file in source_files){
  current_file = readLines(file)
  cat(current_file, sep ="\n")
  cat("\n")
}

sink()