#' Sets the initial state to the steady state
#'
#' The initial a-priori distribution is set to the steady state of
#' the transition matrix.
#'
#' The steady state is computed as the first eigenvector found that
#' the absolute difference of its eigenvalue and one is less than
#' the square root of machine epsilon.
#'
#' The Markov Chain is expected to be irreducible and aperiodic. The
#' first because otherwise the devices would not have freedom of
#' movement. The second because some probabilities from one state
#' to itself are expected to be non zero. This implies that there
#' exists one unique steady state.
#'
#' However, while highly encouraged, it is not enforced and it is
#' possible to define a non irreducible or non aperiodic model. It
#' is not recommendable to use initsteady with such model, as
#' the function only returns one of the steady states.
#'
#'
#' @param x A HMM object.
#'
#' @return The same HMM object of the input with its initial state set
#' to steady state.
#'
#' @examples
#' model <- HMM(2)
#' model <- addtransition(model, c(1,2))
#' model <- addtransition(model, c(2,1))
#' model <- initparams(model)
#' istates(model)
#' model <- initsteady(model)
#' istates(model)
#' (istates(model) %*% getTM(model)) /
#'   sum(istates(model) %*% getTM(model))
#'
initsteady <- function(x) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")

  # Initial steady state calculation
  edecomp <- eigen(t(getTM(x)))
  # The steady state might be not unique!
  i <- which(abs(edecomp$values - 1) < sqrt(.Machine$double.eps))[1]
  states <- Re(edecomp$vectors[,i])
  states <- states / sum(states)

  x$parameters$states <- states

  return(x)
}
