#' This function returns the smooth states

# P : Transition matrix P(will go to i | is in j)
# E : Event observation matrix P(detection event j | is in i)
# e: Sequence of observation events. Since we are taking
# time increase small, it is expected to have mostly missing
# values. The first and last value are expected to have
# observations.
# FS: Filtered/predicted states as returned in fstates

sstates <- function(...) {
  UseMethod("sstates")
}
#' @rdname sstates
sstates.HMM <- function(x, e) {
  fpass <- forward(x,e)
  bpass <- backward(x,e,fpass$scalefactors)
  return(fpass$alpha * bpass)
}
#' @rdname sstates
sstates.matrix <- function(P, E, e, FS) {
  # Allocate the output matrix
  output <- matrix(0, ncol = length(e), nrow = nrow(P))
  # Set initial smoother vector and last smooth state
  svector <- double(ncol(P)) + 1
  output[, length(e)] <- FS[, length(e)]

  # Now run the filter backwards. It is a loop, it should be
  # written in some compiled language for performance.
  for (i in length(e):2) {
    if (!is.na(e[i])) {
      svector <- E[,e[i]] * svector
    }
    svector <- matrix(svector, nrow = 1) %*% P
    svector <- svector / sum(svector)
    output[,i - 1] <- svector * FS[,i - 1]
    output[,i - 1] <- output[,i - 1]/sum(output[,i - 1])
  }
  return(output)
}
