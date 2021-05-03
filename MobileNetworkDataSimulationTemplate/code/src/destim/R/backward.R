#' The backward part of the FB algorithm
#'
#' Calculates the backward probabilities.
#'
#'
#' @return A HMM object.
#'
#' @examples
#' HMM(1L, matrix(c(1L,1L), nrow = 2), EM = matrix(1, nrow = 1))

backward <- function(...) {
  UseMethod("backward")
}
#' @rdname backward
backward.HMM <- function(x,y,sfactors) {
  TM <- getTM(x)
  EM <- emissions(x)
  beta <- matrix(0,nrow = nstates(x), ncol = length(y))
  svector <- rep(1,nstates(x))
  for (i in length(y):1) {
    beta[,i] <- svector
    if (!is.na(y[i]))
      svector <- svector * EM[,y[i]]
    svector <- TM %*% svector
    svector <- svector / sfactors[i]
  }
  return(beta)
}
