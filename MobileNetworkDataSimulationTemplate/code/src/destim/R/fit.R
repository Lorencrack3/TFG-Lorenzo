#' Fits a HMM model
#'
#' Fits the model by maximum likelihood
#'
fit <- function(x, e, sq = FALSE, init = FALSE, method = "constrOptim", ...) {
  if (is.null(x$parameters$reducedparams)) {
    if (is.null(x$parameters$transitions))
      x <- initparams(x)
    x <- minparams(x)
  }
  trmatrix <-
    gettransmatrix(x)[x$parameters$reducedparams$cconstraints, ]
  trmatrix <-
    trmatrix[apply(trmatrix,1,function(trow)
      sum(trow[-length(trow)]^2) != 0), ]
  trmatrix <- funique(trmatrix)
  ui <- trmatrix[, -ncol(trmatrix), drop = FALSE]
  ci <- -trmatrix[, ncol(trmatrix)]
  ui <- rbind(ui, -trmatrix[, -ncol(trmatrix), drop = FALSE])
  ci <- c(ci, trmatrix[, ncol(trmatrix)] - rep(1,length(ci)))
  if (!init)
    ofun <- function(p) {
              rparams(x) <- p
              x <- initsteady(x)
              return(logLik(x,e, sq))
            }
  else
    ofun <- function(p) {
      rparams(x) <- p
      return(logLik(x,e, sq))
    }
  if (method == "constrOptim") {
    if (!require("stats"))
      stop("Package states required for constrOptim method")
    rparams(x) <- stats::constrOptim(rparams(x), ofun, NULL,
                    ui, ci)$par
  }
  else if (method == "donlp2") {
    if (!require("Rdonlp2"))
      stop("Package Rdonlp2 required for donlp2 method")
    rparams(x) <-
      Rdonlp2::donlp2NLP(rparams(x), ofun,
        par.lower = rep(0, length(rparams(x))),
        par.upper = rep(1, length(rparams(x))),
        ineqA = trmatrix[, -ncol(trmatrix)],
        ineqA.lower = -trmatrix[, ncol(trmatrix)],
        ineqA.upper = 1 - trmatrix[, ncol(trmatrix)], ...)$solution
  }
  else if (method == "solnp") {
    if (!require("Rsolnp"))
      stop("Package Rsolnp required for solnp method")
    rparams(x) <- Rsolnp::solnp(rparams(x), ofun,
                  ineqfun = function(MAT) trmatrix %*% c(MAT, 1),
                  ineqLB = rep(0, nrow(trmatrix)),
                  ineqUB = rep(1, nrow(trmatrix)), ...)$pars
  }
  else stop(paste0("Method unknown: ", method))
if (!init)
  x <- initsteady(x)
return(x)
}
