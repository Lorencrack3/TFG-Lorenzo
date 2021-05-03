#' @title Generates random values from a Poisson multinomial distribution.
#'
#' @description Generates random values from a Poisson multinomial distribution.
#'   This function is only called internally. It cannot be called directly by
#'   users.
#'
#' @param n The number of random values to generate.
#'
#' @param prob.dt A data.table object with the following columns:\code{device,
#'   cell, devCount, prob}.
#'
#' @return A matrix object with the random values generated according to a
#'   Poisson multinomial distribution.
#'
#'
#' @import data.table
#' @import extraDistr
rNnet_Event <- function(n, prob.dt){
  
  
  if (!all(c('device', 'region', 'devCount', 'prob') %in% names(prob.dt))) {
    
    stop('[rNnet_Event] prob.dt must have columns device, cell, devCount, prob.\n')
  }
  
  probSums <- prob.dt[, list(totalProb = sum(prob)), by = 'device']$totalProb
  if (!all(abs(probSums - 1) < 1e-8)) {
    #stop('[rNnet_Event] The sum of probabilities per device is not 1.\n')
  }
  
  x1 <- prob.dt[
      , c('device', 'region', 'devCount', 'prob'), with = FALSE][
        , categories := paste0(region, '-', devCount)]
  x2 <- x1[
    , list(category = rcat(n, prob, categories)), by = 'device'][
      , nSim := 1:n, by = 'device']
  x3 <- dcast(x2, device ~ nSim, value.var = 'category')
  cells <- names(nIndividuals(x3[[2]]))
  x4 <- t(as.matrix(x3[
    , lapply(.SD, nIndividuals), .SDcols = names(x3)[-1]]))
  dimnames(x4)[[2]] <- cells
  return(x4)
}
