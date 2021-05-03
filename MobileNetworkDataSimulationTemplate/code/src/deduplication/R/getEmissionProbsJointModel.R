#' @title Builds the emission probabilities for the joint HMM.
#'
#'
#' @description Builds the emissions probabilities needed for the joint HMM used to estimate the posterior location
#'   probabilitities.
#'
#' @param emissionProbs the emission probabilities (the location probabilities) computed by calling
#'   \code{getEmissionProbs()} for each individual device.
#'
#'
#' @return Returns a matrix with the joint emission probabilities for the HMM. The number of rows equals the number
#'   tiles and the number of columns equals the number of combinations between antenna IDs. Before the combination
#'   between antenna IDs are build, the NA value is added to the list of antenna IDs. An element in this matrix
#'   represents the transition probability from an antenna to another, computed for each tile in the grid.
#'
#'
#' @import Matrix
#'
#' @export
getEmissionProbsJointModel <- function(emissionProbs) {
  
  emissionProbs <-
    cbind(emissionProbs, 'NA' = rep(1, nrow(emissionProbs)))
  jointEmissionProbs <- NULL
  for (j in 1:ncol(emissionProbs)) {
    A <- emissionProbs[, j] * emissionProbs
    colnames(A) <-
      paste0(colnames(emissionProbs)[j], "-", colnames(emissionProbs))
    jointEmissionProbs <- cbind(jointEmissionProbs, A)
  }
  colToRem <- ncol(jointEmissionProbs)
  jointEmissionProbs <- jointEmissionProbs[,-colToRem]
  
  return(jointEmissionProbs)
  
}
