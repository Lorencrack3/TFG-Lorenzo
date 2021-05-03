#' Functions to modify a HMM object
#'
`rparams<-` <- function(x, value) {
  if (is.null(x$parameters$reducedparams)) {
    if (is.null(x$parameters$states))
      x <- initparams(x)
    x <- minparams(x)
  }
  x$parameters$reducedparams$params <- value
  x$parameters$transitions <-
    (x$parameters$reducedparams$transmatrix %*%
      matrix(c(x$parameters$reducedparams$params, 1), ncol = 1)) [,]
  return(x)
}

`setsnames<-` <- function(x,value) {
  if (length(value) != nstates(x))
    stop("Wrong number of states")
  x$states$names <- value
  return(x)
}

`emissions<-` <- function(x, value) {
  if (nrow(value) != nstates(x))
    stop("The number of rows does not match the number of states.")
  x$emissions <- value
  return(x)
}

`istates<-` <- function(x,value) {
  if (length(value) != nstates(x))
    stop("Wrong number of states")
  x$parameters$states <- value
  return(x)
}
