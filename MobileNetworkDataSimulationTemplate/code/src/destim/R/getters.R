#' Number of states getter
#'
nstates <- function(x) {
  UseMethod("nstates")
}
#' @rdname nstates
nstates.HMM <- function(x) return(length(x$states$names))

#' Number of transition parameters getter
#'
ntransitions <- function(x) {
  UseMethod("ntransitions")
}
#' @rdname ntransitions
ntransitions.HMM <- function(x) return(ncol(x$transitions))
#' Number of constraints getter
#'
nconstraints <- function(x) {
  UseMethod("nconstraints")
}
#' @rdname nconstraints
nconstraints.HMM <- function(x) return(nrow(x$constraints))
#' Matrix of constraints getter
#'
constraints <- function(x) {
  UseMethod("constraints")
}
#' @rdname constraints
constraints.HMM <- function(x) return(x$constraints)
#' Probability of transitions parameters getter
#'
ptransition <- function(x) {
  UseMethod("ptransition")
}
#' @rdname ptransition
ptransition.HMM <- function(x)
  return(as.numeric(x$parameters$transitions))
#' Initial state parameters getter
#'
istates <- function(x) {
  UseMethod("istates")
}
#' @rdname istates
istates.HMM <- function(x) return(x$parameters$states)
#' Transitions getter
#'
transitions <- function(x) {
  UseMethod("transitions")
}
#' @rdname transitions
transitions.HMM <- function(x) return(x$transitions)
#' Reduced parameters getter
rparams <- function(x) {
  UseMethod("rparams")
}
#' @rdname rparams
rparams.HMM <- function(x) return(x$parameters$reducedparams$params)
#' Transformation matrix getter
gettransmatrix <- function(x) {
  UseMethod("gettransmatrix")
}
#' @rdname gettransmatrix
gettransmatrix.HMM <- function(x)
  return(x$parameters$reducedparams$transmatrix)
#' Emissions matrix getter
#'
emissions <- function(x) {
  UseMethod("emissions")
}
#' @rdname emissions
emissions.HMM <- function(x) return(x$emissions)
#' Transitions matrix getter.
#'
#' Returns the transitions matrix from a HMM.
#'
#' The HMM object contains three fields. States contains information
#' about the states. Transitions contains a list of matrices with
#' two rows. Each column represents a transition with non-zero
#' probability. Transitions in the same matrix have same probability.
#' Emissions is a matrix with the emission probabilities. These are
#' considered fixed.
#'
#' @param x  An HMM object.
#'
#' @return The transition matrix of the HMM.
#'
#' @examples
#' HMM(1L, list(matrix(c(1L,1L), nrow = 2)), EM = matrix(1, nrow = 1))

getTM <- function(x) {
  UseMethod("getTM")
}
#' @rdname getTM
getTM.HMM <- function(x) {
  if (is.null(x$parameters$transitions))
     stop(paste0("[destim::getTM] Parameters of the model are ",
      "required to get the transitions matrix"))
  TM <- matrix(0,nrow = nstates(x), ncol = nstates(x))
  for (i in 1:ntransitions(x))
      TM[x$transitions[1L,i],
          x$transitions[2L,i]] <- x$parameters$transitions[i]
  return(TM)
}
