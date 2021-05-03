#' @title Computes tau coefficient
#'
#' @description Computes tau, a coefficient needed for the calculation of the population counts at time t>t0.
#'
#' @param nnetODFileName the name of the file where the population moving from one region to another is stored. This is
#'   an output of the \code{aggregation} package. Since this file could be large, it is stored compressed in zip format.
#'
#' @param zip If TRUE the file where where the population moving from one region to another is stored is a zipped csv
#'   file, otherwise it is  simple csv file.
#'
#' @import data.table
#' @import utils
#' @export
computeTau <- function(nnetODFileName, zip) {
    
    if (!file.exists(nnetODFileName))
        stop(paste0(nnetODFileName, " does not exists!"))
    
    if(zip ==TRUE) {
      nnetUZ<-unzip(nnetODFileName, exdir=tempdir())
    } else { 
      nnetUZ = nnetODFileName
    }
    nNetOD <- fread(
        nnetUZ,
        sep = ',',
        header = TRUE,
        stringsAsFactors = FALSE
    )
    
    NnetOD_totalFrom <- nNetOD[
        , list(Nnet_totalFrom = sum(Nnet)), by = c('time_from', 'time_to', 'region_from', 'iter')]
    
    
    tauOD <- merge(nNetOD, NnetOD_totalFrom, 
      by = c('time_from', 'time_to', 'region_from', 'iter'))[
      , tau_Nnet := Nnet / Nnet_totalFrom][
        Nnet_totalFrom == 0, tau_Nnet := 0][
      , Nnet_totalFrom := NULL][
      , Nnet := NULL]
    
    return (tauOD)
}