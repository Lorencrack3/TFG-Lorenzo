#' Dissimilarity
#'
#' This is an approximate measure
#'
#'

dissim <- function(...) {
  UseMethod("dissim")
}

#' @rdname dissim
dissim.HMM <- function(x, y1, y2) {
  E <- sapply(y2, function(z)
    if (!is.na(z)) emissions(x)[, z]
    else rep(1, nstates(x)))
  output <- apply(E * sstates(x, y1), 2, sum)
  return(-sum(log(output)) / sum(!is.na(y2)))
}
