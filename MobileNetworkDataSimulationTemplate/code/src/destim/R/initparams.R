#' A simple Initializer for HMM objects
#'
#' Set initial parameters for an HMM object, as specified.
#'
#' The field parameters of the HMM object, which includes both
#' initial state and transition probabilities, is initialized at
#' random.
#'
#' The initial states probabilities are set to an uniform
#' (0,1) distribution and then divided by their sum.
#'
#' The initial probabilities of transition are also set to an uniform
#' (0,1) and in this case, projected on the constrained space. After
#' the projection some probability might result greater than one or
#' less than zero. Those probabilities are then set to uniform (0,1)
#' again and the process is repeated until all probabilities of
#' transition are in (0,1) and the constraints are satisfied.
#'
#' @param x A HMM object.
#'
#' @return An initialized HMM object.
#'
#' @examples
#' model <- HMMrectangle(2,2)
#' model <- initparams(model)
#' range(constraints(model) %*% c(ptransition(model), -1))
#'
initparams <- function(x) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")

  # Sets random parameters before enforcing constraints
  states <- runif(nstates(x))
  transitions <- runif(ntransitions(x))

  # Enforce constraints: we calculate the closest point in the
  # subspace using Lagrange multipliers.
  # It is expected that we have many constraints that are
  # equalities between transition probabilities. Those constraints
  # allow us to remove one of them from the equations.
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
  # Now we only need non equality constraints in the system
  MCT <- CT[!eqcon, , drop = FALSE]
  MCT <- cbind(MCT[, -ncol(MCT), drop = FALSE] %*% trmatrix,
               MCT[,ncol(MCT), drop = FALSE])
  MCT <- funique(MCT) # Just in case
  # Now we get the system of equations to obtain the closest
  # point in the constrained space.
  eqsys <- cbind(matrix(0,nrow = nrow(MCT),
                        ncol = nrow(MCT)),
                 MCT)
  eqsys <- rbind(eqsys,
                 cbind(t(MCT[, -ncol(MCT), drop = FALSE]),
                     diag(ncol(trmatrix)) * 2,
                     matrix(2 * transitions[cconstraints], ncol = 1)))
  while (TRUE) {
  transitions <-
    solve(eqsys[,1:nrow(eqsys)], eqsys[,ncol(eqsys)])
  transitions <- transitions[-1:(-nrow(MCT))]
  # If all transition probabilities are between 0 and 1 now,
  # we are done.
  if (all(transitions >= 0 & transitions <= 1))
    break
  # Otherwise, set again random values for those and repeat.
  transitions[transitions < 0] <- runif(sum(transitions < 0))
  transitions[transitions > 1] <- runif(sum(transitions > 1))
  eqsys[-1:(-nrow(MCT)), ncol(eqsys)] <- 2 * transitions
  }
  states <- states / sum(states)
  x$parameters <- list(states = states,
                       transitions = trmatrix %*% transitions)
  return(x)
}
