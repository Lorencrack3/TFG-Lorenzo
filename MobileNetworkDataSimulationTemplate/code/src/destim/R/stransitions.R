#' This function returns Xi
#'
#'

# P : Transition matrix P(will go to i | is in j)
# E : Event observation matrix P(detection event j | is in i)
# e: Sequence of observation events. Since we are taking
# time increase small, it is expected to have mostly missing
# values. The first and last value are expected to have
# observations.
# FS: Filtered/predicted states as returned in fstates


scpstates <- function(...) {
  UseMethod("scpstates")
}
#' @rdname scpstates
scpstates.HMM <- function(x, e) {
  output <- matrix(0, nrow = nstates(x) ** 2, ncol = length(e) - 1)
  TM <- getTM(x)
  EM <- emissions(x)
  fpass <- forward(x,e)
  bpass <- backward(x,e,fpass$scalefactors)
  for (i in 2:length(e)) {
    if (is.na(e[i]))
      output[, i - 1] <-
        (TM * (matrix(fpass$alpha[, i - 1], ncol = 1) %*%
        matrix(bpass[, i], nrow = 1)) )[1:nrow(output)] /
        fpass$scalefactors[i]
    else
      output[, i - 1] <-
        (TM * (matrix(fpass$alpha[, i - 1], ncol = 1) %*%
        matrix(bpass[, i] * EM[,e[i]], nrow = 1)) )[1:nrow(output)] /
        fpass$scalefactors[i]
  }
  return(output)
}
