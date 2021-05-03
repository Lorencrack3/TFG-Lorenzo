#' @title Returns the mode of delta distribution.
#'
#' @description Returns the mode of the deltaX or deltaY distribution.
#'
#' @param deltaDistribution A data.table object that could be Delta X or Delta Y distribution. The table has two columns:
#'   delta and p. It is obtainded from calling \code{buildDeltaProb()} function.
#'
#' @return Rhe mode of the delta distribution.
#'
#' @import data.table
#'
#' @export
modeDelta <-function(deltaDistribution) {
  
  return(deltaDistribution[which.max(deltaDistribution$p), delta])
  
}
