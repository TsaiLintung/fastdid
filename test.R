f <- function() {
  pause(1)
  g()
  h()
}
g <- function() {
  pause(1)
  h()
}
h <- function() {
  pause(2)
}