build_source <- function(ver){
  source_files <- list.files("R", include.dirs = FALSE, full.names = TRUE)

  dep <- c("data.table", "stringr", "BMisc", "collapse", "dreamerr", "parglm")
  
  sink(paste0("development/fastdid_", str_replace_all(ver, "\\.", "_"), ".R"))
  cat(paste0("#", as.character(Sys.Date()), "\n"))
  cat(paste0("message('loading fastdid source ver. ver: ", ver, " date: " , as.character(Sys.Date()), "')\n"))
  cat(paste0("require(", dep, ");" |> str_flatten(), "\n"))
  for(file in source_files){
    current_file = readLines(file, warn = FALSE)
    cat(current_file, sep ="\n")
    cat("\n")
  }
  
  sink()
}
