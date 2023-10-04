


set_max_thread <- function(){
  setDTthreads(0)
  options(kit.nThread = getDTthreads())
  setFixest_nthreads(getDTthreads())
}

reverse_col <- function(x){
  return(x[,ncol(x):1])
}