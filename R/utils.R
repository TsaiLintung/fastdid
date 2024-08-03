set_max_thread <- function(){
  data.table::setDTthreads(0)
  options(kit.nThread = getDTthreads())
}

reverse_col <- function(x){
  return(x[,ncol(x):1])
}

#g11 <- c(1,1,2,3,4)
#g22 <- c(2,1,3,2,4)
#GG <- as.factor(paste0(g11, ".", g22))

g1 <- function(GG){
  if(is.numeric(GG)){return(GG)}
  return(as.numeric(str_split_i(as.character(GG), "-", 1)))
}
g2 <- function(GG){
  return(as.numeric(str_split_i(as.character(GG), "-", 2)))
}
ming <- function(GG){
  if(is.numeric(GG)){return(GG)}
  else {pmin(g1(GG), g2(GG))}
}