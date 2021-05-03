#' Reparametrizes a HMM model with a minimal set of parameters.
#'
#' Finds a minimal set of parameters that fully determine the
#' probabilities of transition.
#'
#' This function avoids to solve a high dimensional optimization
#' problem with many constraints, parametrizing the probabilities
#' of transition with as few parameters as possible: the number of
#' degrees of freedom.
#'
#' A pivoted QR decomposition of the constraints matrix is done,
#' to get both the free parameters and the matrix that transforms
#' them back into the probabilities of transition.
#'
#' Many constraints are expected to be equalities between two
#' probabilities of transition, so the function is optimized for
#' this special kind of constraints.
#'
#' @param x A HMM object.
#'
#' @return The same HMM object of the input with some additional
#' fields that store the new parameters and the matrix that
#' transforms these parameters in the probabilities of transition.
#'
#' @seealso \link{rparams}, \link{rparams<-}
#'
#' @examples
#' model <- HMMrectangle(2,2)
#' model <- initparams(model)
#' ntransitions(model)
#' nconstraints(model)
#' model <- minparams(model)
#' rparams(model)
#' range(ptransition(model) -
#'   gettransmatrix(model) %*% c(rparams(model), 1))
#'
minparams <- function(x) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")
  trmatrix <- matrix(as.logical(diag(ntransitions(x))),
                     ncol = ntransitions(x))
  CT <- constraints(x)
  eqcon <- apply(CT, 1, function(irow)
    (irow[ntransitions(x) + 1] == 0) &&
      (sum(irow != 0) == 2) &&
      (sum(irow) == 0))
  EQCT <- CT[eqcon, 1:ntransitions(x), drop = FALSE]
  trmatrixcols <- ncol(trmatrix)
  for (i in seq_len(nrow(EQCT))) {
    variables <- match(c(1, -1), EQCT[i, ], incomparables = 0)
    columns <- c(match(TRUE, trmatrix[variables[1], ]),
                 match(TRUE, trmatrix[variables[2], ]))
    columns <- sort(columns)
    trmatrix[,columns[1]] <- trmatrix[,columns[1]] |
                               trmatrix[,columns[2]]
    trmatrix[, columns[2]] <-
     trmatrix[, trmatrixcols, drop = FALSE]
    trmatrixcols <- trmatrixcols - 1
  }
  trmatrix <- trmatrix[,1:trmatrixcols]
  cconstraints <-
    apply(trmatrix, 2, function(v) which(v)[1])
  trmatrix <- matrix(as.numeric(trmatrix), ncol = ncol(trmatrix))
  MCT <- CT[!eqcon, , drop = FALSE]
  MCT <- cbind(MCT[, -ncol(MCT), drop = FALSE] %*% trmatrix,
               MCT[,ncol(MCT), drop = FALSE])
  qrdecomp <- qr(MCT[, -ncol(MCT), drop = FALSE])
  MCT <- cbind(qr.R(qrdecomp), -t(qr.Q(qrdecomp)) %*% MCT[, ncol(MCT)])
  MCT <- -solve(MCT[1:qrdecomp$rank,1:qrdecomp$rank],
                MCT[1:qrdecomp$rank,-(1:qrdecomp$rank),
                    drop = FALSE])
  MCT <- rbind(MCT, cbind(diag(length(qrdecomp$pivot) - nrow(MCT)),
                          rep(0,length(qrdecomp$pivot) - nrow(MCT))))
  trmatrix <- trmatrix[, qrdecomp$pivot] %*% MCT
  if (!is.null(x$parameters$reducedparams))
    x$parameters <- x$parameters[c("states", "transitions")]
  x$parameters <- c(x$parameters,list(reducedparams =
                         list(transmatrix = trmatrix,
                              cconstraints = cconstraints,
                              params = NULL)))
  if (!is.null(x$parameters$transitions)) {
    x$parameters$reducedparams$params <- x$parameters$transitions[
      sapply(1:(ncol(trmatrix) - 1), function(x)
        which(apply(trmatrix ** 2, 1, function(y) sum(y) == y[x]))[1])
      ]
    x$parameters$transitions <- (trmatrix %*%
      matrix(c(x$parameters$reducedparams$params, 1), ncol = 1)) [,]
  }
  return(x)
  }
