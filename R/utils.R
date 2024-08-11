

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