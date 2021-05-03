#' @title Generates random values from a Poisson multinomial distribution needed
#'   for the origin destination matrices.
#'
#' @description Generates random values from a Poisson multinomial distribution
#'   needed for the origin destination matrices. This function is only called
#'   internally and it cannot by called directly by users.
#'
#' @param n The number of random values to generate.
#'
#' @param prob.dt A data.table object with the following columns:\code{device,
#'   cell, devCount, prob}.
#'
#' @return A matrix object with the random values generated according to a
#'   Poisson multinomial distribution.
#'
#' @import data.table
#' @import extraDistr
#' @include nIndividuals3.R
#' @export
rNnetJoint_Event <- function(n, prob.dt, cellNames){
  
  
  if (!all(c('device', 'region_from', 'region_to', 'devCount', 'prob') %in% names(prob.dt))) {
    
    stop('[rNnetCond_E] prob.dt must have columns device, region_from, region_to, devCount, prob.\n')
  }
  
  #probSums <- prob.dt[, list(totalProb = sum(prob)), by = c('device', 'region_from', 'region_to')]$totalProb
  #if (!all(abs(probSums - 1) < 1e-8)) {
  #  
  #  stop('[rNnetCond_E] The sum of probabilities per device is not 1.\n')
  #}
  
  x1 <- prob.dt[
    , c('device', 'region_from', 'region_to', 'devCount', 'prob'), with = FALSE][
      , categories := paste0(region_from, '_', region_to, '-', devCount)]

  x2 <- x1[
    , list(category = rcat(n, prob, categories)), by = c('device')][
    , nSim := 1:n, by = c('device')][
    , c('region_fromto') := tstrsplit(category, split = '-', keep = 1)]
  x3 <- dcast(x2, device + region_fromto ~ nSim, value.var = 'category')
  x3 <- x3[
    , lapply(.SD, function(x){x[!is.na(x)]}), by = 'device', .SDcols = names(x3)[-c(1:2)]]
  tempPairs <- data.table(expand.grid(cellNames, cellNames))[
    , cellPairs := do.call(paste, c(.SD, sep = "_"))]
  cellPairs <- tempPairs$cellPairs
  x4 <- x3[
    , lapply(.SD, nIndividuals3, cellNames = cellNames), .SDcols = names(x3)[-1]][
    , cellPairs := cellPairs][
    , c('region_from', 'region_to') := tstrsplit(cellPairs, split = '_')][
    , cellPairs := NULL]
  
  x4.molten <- melt(x4, id.vars = c('region_from', 'region_to'),
                    value.factor = FALSE, value.name = 'Nnet',
                    variable.name = 'iter')
  setcolorder(x4.molten, c('region_from', 'region_to', 'Nnet', 'iter'))
  return(x4.molten)  
  
}

