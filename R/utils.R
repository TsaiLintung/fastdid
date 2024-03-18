set_max_thread <- function(){
  data.table::setDTthreads(0)
  options(kit.nThread = getDTthreads())
}

reverse_col <- function(x){
  return(x[,ncol(x):1])
}


release <- function(p){
  for(name in names(p)){
    assign(name, p[[name]], envir = parent.frame()) #assign it to the envir where the func is called
  }
}

gather <- function(...){
  names <- substitute(list(...)) |> deparse() |>
    str_remove_all("\\)|list\\(") |> str_flatten() |> str_split_1(",") |> str_trim()
  values <- list(...)
  p <- list()
  
  for(i in 1:length(names)){
    #p[[names[[i]]]] <- values[[i]]
    p <- c(p, list(values[[i]])) #for nulls
  }
  names(p) <- names
  return(p)
}
