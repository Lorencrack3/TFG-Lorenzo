#' Class constructor for Hidden Markov models
#'
#' Creates an HMM object, as specified.
#'
#' The HMM object contains five fields: states, transitions,
#' constraints, emisions and parameters.
#'
#' The field states contains a character vector with the names of
#' the states. If the constructor is given a number S, it sets the
#' names as follows: as.character(1:S). An additional field, called
#' coordinates is provided too, were in the future the geolocation of
#' the states will be specified. Note that state determines geolocation,
#' but different states might share geolocation.
#'
#' The field transitions contain a matrix which is a list of the
#' transitions with non-zero probability. It is a two row integer matrix
#' where each column represents the transition from first row state
#' to second row state. The states are referenced in the same order as
#' they appear in field states. While (number of states)^2 transitions
#' are possible, a much smaller number is expected. It defaults to still
#' transitions for all states.
#'
#' The field constraints is the augmented matrix of the system of
#' linear equalities that the model must fulfill. The variables of the
#' system correspond to the probabilities of transition, in the same
#' order as in field transitions. The constraints that restrict the
#' transition matrix to be stochastic are treated automatically.
#'
#' The field emissions consists in a matrix that contains the emission
#' probabilities, where the number of rows is the number of states and
#' each column correspond to a possible output. Unlike usual, the
#' emission probabilities are fixed.
#'
#' The field parameters contain additional information about the
#' probabilities of transition and the initial state of the model.
#' Also some auxiliar information to reduce the number of parameters
#' of the model. See initparams, minparams and initsteady.
#'
#'
#'
#' @param S  Number or names of states. It can be either a numeric or a
#' character.
#' @param TL Matrix of integers that lists non-zero transitions.
#' The matrix corresponds to the field transitions of the object
#' (see details).
#' @param CT Matrix of constraints. It corresponds to the field
#' constraints of the object (see details).
#' @param EM Matrix of emissions. It corresponds to the field
#' emissions of the object (see details).
#'
#' @return A HMM model object.
#'
#' @seealso \link{initparams}, \link{minparams}, \link{initsteady}
#'
#' @examples
#' model1 <- HMM(5)
#' model2 <- HMM(c("a","b","c"),
#'             TL = matrix(c(1, 1,
#'                           1, 2,
#'                           2, 1,
#'                           2, 2,
#'                           2, 3,
#'                           3, 2,
#'                           3, 3), nrow = 2))
#' nstates(model1)
#' ntransitions(model1)
#' nstates(model2)
#' ntransitions(model2)

HMM <- function(...) {
  UseMethod("HMM")
}
#' @rdname HMM
HMM.integer <- function(S, TL, CT, EM = NULL) {

  if (missing(TL))
    TL = matrix(c(1:S,1:S), nrow = 2, byrow = TRUE)
  if (nrow(TL) != 2)
    stop("The list of transitions must contain exactly two rows.")
  if ((min(TL) < 1) || max(TL) > S)
    stop("The states are referenced by numbers from 1 to S")
  # Remove possible duplicates
  TL <- t(unique(t(TL)))

  BCT <- t(sapply(1:S, function(x) c(TL[1,] == x, 1.0)))
  if (missing(CT))
    CT <- BCT
  else {
    if (ncol(CT) != ncol(TL) + 1)
      stop(paste0("The number of columns of the constraints matrix ",
                   "does not match the number of transitions."))
    CT <- rbind(CT, BCT)
  }
  # Remove possible duplicates
  # TODO

  output <- list(states = list(names = as.character(1:S),
                               coordinates = NULL),
                 transitions = TL, constraints = CT,
                 emissions = EM, parameters = list(transitions = NULL,
                                                   states = NULL))
  attr(output,"class") <- "HMM"
  return(output)
}
#' @rdname HMM
HMM.numeric <- function(S, ...) return(HMM(as.integer(S), ...))
#' @rdname HMM
HMM.character <- function(S, ...) {
  output <- HMM(length(S), ...)
  setsnames(output) <- S
  return(output)
}
