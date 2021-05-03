#' Adds a transition to the model
#'
#' The specified transition is added to the model as a
#' transition with non zero probability.
#'
#' Since the transition probabilities from the initial state of the
#' newly specified transition still have to sum up to one, the
#' correspondent constraint is modified accordingly.
#'
#' @param x A HMM object.
#' @param t The transition, as a two dimensional integer vector.
#' The first element is the number of the initial state and the
#' second one the number of the final state.
#'
#' @return A HMM object similar to the input but with the additional
#' transition.
#'
#' @seealso \link{HMM}, \link{addconstraint}
#'
#' @examples
#' model <- HMM(3)
#' model <- addtransition(model, c(1,2))
#' model <- addtransition(model, c(2,3))
#' model <- addtransition(model, c(3,1))
#' transitions(model)
#' constraints(model)
#'
addtransition <- function(x,t) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")
  if (length(t) < 2)
    stop("Initial and final state must be specified.")
  if (length(t) > 2) {
    warning(paste0("Only initial and final state are required. ",
                   "Remaining elements of the vector are ignored."))
    t <- t[1:2]
  }
  if (any(apply(transitions(x), 2, function(z) all(z == t))))
    stop("Transition already stated as having non zero probability.")
  TL <- x[["transitions"]]
  CT <- x[["constraints"]]
  CTrow <- which(sapply(1:nrow(CT), function(z) {
                return(all(as.numeric(TL[1,] == t[1]) ==
                             CT[z, -ncol(CT)]))
    }))
  TL <- cbind(matrix(t, ncol = 1), TL)
  CT <- cbind(matrix(0,nrow = nrow(CT), ncol = 1), CT)
  CT[CTrow, 1] <- 1
  x[["transitions"]] <- TL
  x[["constraints"]] <- CT
  x[["parameters"]] <- list(transitions = NULL, states = NULL)
  return(x)
}
