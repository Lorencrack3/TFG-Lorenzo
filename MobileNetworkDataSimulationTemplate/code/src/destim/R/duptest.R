#' Likelihood ratio test
#'
#' The likelihood ratio test is used to check whether two sets of observed events
#' come from two devices carried by the same person.
#'
#'

duptest <- function(x, y1, y2) {
  UseMethod("duptest")
}

#' @rdname duptest
duptest.HMM <- function(x, y1, y2) {
  if (length(y1) != length(y2))
    stop("Both sets of observed events have to ")
  ll <- 0
  set1 <- numeric(2 * length(y1))
  set1[seq(2, length(set1), 2)] <- NA_integer_
  set1[seq(1, length(set1), 2)] <- y1
  set2 <- numeric(2 * length(y1))
  set2[seq(1, length(set1), 2)] <- NA_integer_
  set2[seq(2, length(set1), 2)] <- y2
  x <- initparams(x)
  x <- minparams(x)
  x <- fit(x, set1)
  ll <- ll - logLik(x, set1)
  x <- initparams(x)
  x <- minparams(x)
  x <- fit(x, set2)
  ll <- ll - logLik(x, set2)
  set1[seq(2, length(set1), 2)] <- y2
  x <- initparams(x)
  x <- minparams(x)
  x <- fit(x, set1, sq = TRUE)
  ll <- ll + logLik(x, set1, sq = TRUE)
  ll <- 2 * ll
  return(list(LR = ll, df = length(set1) + 1 + length(rparams(x))))
}
