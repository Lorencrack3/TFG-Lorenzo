#' @title Computes the number of individuals moving from a region to another as
#'   a sum of multinomial variables.
#'
#' @description Computes the number of individuals moving from a region to
#'   another as a sum of multinomial varibales. It is an internal function that
#'   is called during the process of generating the Poisson multinomial
#'   distributed random values. This function cannot be directly called by
#'   users.
#'
#' @param categories a list of categories in the form \code{region-devCount}.
#'
#' @param cellNames A vector with the names of the regions.
#'
#' @return The sum of the multinomial variates.
#'
#' @import data.table
nIndividuals3 <- function(categories, cellNames){
  
  catPerDevice <- Reduce(rbind, lapply(as.character(categories), function(str) {
    
    strsplit(str, split = '-')[[1]]
  }))
  nDev <- length(categories)
  cells <- sort(unique(catPerDevice[, 1]))
  tempPairs <- data.table(expand.grid(cellNames, cellNames))[
    , cellPairs := do.call(paste, c(.SD, sep = "_"))]
  cellPairs <- tempPairs$cellPairs
  if (!all(cells %in% cellPairs)) stop('[nIndividuals] cells')
  output <- numeric(length(cellPairs))
  names(output) <- cellPairs
  for (dev_index in 1:nDev){
    
    cell <- catPerDevice[dev_index, 1]
    devCount <- as.numeric(catPerDevice[dev_index, 2])
    tempVec <- numeric(length(cellPairs))
    names(tempVec) <- cellPairs
    tempVec[as.character(cell)] <- devCount
    output <- output + tempVec    
  }
  
  return(output)
}
