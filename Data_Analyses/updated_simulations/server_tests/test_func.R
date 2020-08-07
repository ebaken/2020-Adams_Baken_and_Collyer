test_func <- function(x,y) {
  a <- rnorm(x)
  b <- rnorm(y)
  length(a) - max(b)
}