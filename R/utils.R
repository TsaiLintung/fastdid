set_max_thread <- function(){
  data.table::setDTthreads(0)
  options(kit.nThread = getDTthreads())
}

reverse_col <- function(x){
  return(x[,ncol(x):1])
}

#locked list
#from: https://stackoverflow.com/questions/58776481/make-r-function-return-a-locked-immutable-list
#' @export
`[[<-.locked` <- function(value) {stop("Can't assign into locked object")}
#' @export
`[<-.locked` <- function(value) {stop("Can't assign into locked object")}
#' @export
`$<-.locked` <- function(value) {stop("Can't assign into locked object")}
# a <- list(b = 1, c = 2)
# class(a) <- c("locked")