#' @title Builds the generic HMM model.
#'
#' @description Builds the generic HMM model using the emission probabilities given by \code{getEmissionProbs()}.
#'
#' @param nrows Number of rows in the grid.
#'
#' @param ncols Number of columns in the grid.
#'
#' @param emissionProbs A matrix with the event location probabilities. The number of rows equals the number of tiles in
#'   the grid and the number of columns equals the number of antennas. This matrix is obtained by calling
#'   \code{getEmissionProbs()} function.
#'
#' @param initSteady If TRUE the initial apriori distribution is set to the steady state of the transition matrix, if
#'   FALSE the apriori distribution should be given as a parameter.
#'
#' @param aprioriProb The apriori distribution for the HMM model. It is needed only if initSteady is FALSE.
#'
#' @return Returns an HMM model with the initial apriori distribution set to the steady state of the transition matrix
#'   or to the value given by \code{aprioriProb} parameter.
#'
#' @import destim
#'
#' @export
getGenericModel <-  function(nrows, ncols, emissionProbs, initSteady = TRUE, aprioriProb = NULL) {
  model <- HMMrectangle(nrows, ncols)
  emissions(model) <- emissionProbs
  model <- initparams(model)
  model <- minparams(model)
  
  if( initSteady && !is.null(aprioriProb)) {
    stop("getGenericModel: either initSteady is TRUE and aprioriProb is NULL or initSteady is FALSE and aprioriProb in 
         not NULL")
  }
  
  if( !initSteady && is.null(aprioriProb)) {
    stop("getGenericModel: either initSteady is TRUE and aprioriProb is NULL or initSteady is FALSE and aprioriProb in 
         not NULL")
  }
  
  if (initSteady) {
    model <- initsteady(model)
  } else {
    if(is.null(aprioriProb) ) {
      stop("getGenericModel: if initSteady is FALSE then you should specify the apriori probability for the HMM model!")
    }
    else {
      if(abs(sum(aprioriProb) - 1) > .Machine$double.eps)
        stop("getGenericModel: aprioriProb should sum up to 1!")
      else
        istates(model)<-aprioriProb
    }
  }
  return (model)
}