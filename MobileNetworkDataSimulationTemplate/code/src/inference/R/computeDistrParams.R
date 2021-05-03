#' @title Computes the parameters of the population counts distributions.
#'
#'
#'
#' @param omega
#'
#' @param popRegFileName The name of the file with the population counts for each region taken from a population
#'   register. It has 2 columns: \code{region, Nreg}.
#'
#' @param pntRateFileName The name of the file with the penetration rates for each region. It has 2 columns:
#'   \code{region, pntRate}.
#'
#' @param regsFileName The name of the .csv file defining the regions. It has two columns: \code{ tile, region}. The
#'   first column contains the IDs of each tile in the grid while the second contains the number of a region. This file
#'   is defined by the user and it can be created with any text editor. It is required only for the state process
#'   negative binomial distribution.
#'
#' @param gridFileName The name of the .csv file with the grid parameters.  It is required only for the state process
#'   negative binomial distribution.
#'
#' @param re_bias The value of the relative bias for the population density of each region. The default value is 0.
#' 
#' @param cv The coefficient of variation for the population density of each region. The default value is 0.
#' 
#' @return A data.table object with the following columns  \code{region, omega1, omega2, pnrRate, regionArea_km2, Nreg,
#'   dedupPntRate, alpha, beta}. If regsFileName and gridFileName are not NULL the result will have 3 more
#'   columns:\code{region, omega1, omega2, pnrRate, regionArea_km2, Nreg, dedupPntRate, alpha, beta, theta, zeta, Q}. They
#'   are needed only for the state process negative binomial distribution.
#'
#'
#' @import data.table
#' @export
computeDistrParams <- function(omega, popRegFileName, pntRateFileName, regsFileName = NULL, gridFileName = NULL, rel_bias = 0, cv = 1e-5) {
    
    if (!file.exists(popRegFileName))
        stop(paste0(popRegFileName, " does not exists!"))
    
    pop_reg <- fread(
        popRegFileName,
        sep = ',',
        header = TRUE,
        stringsAsFactors = FALSE
    )
    
    if (!file.exists(pntRateFileName))
        stop(paste0(pntRateFileName, " does not exists!"))
    
    pntRate <- fread(
        pntRateFileName,
        sep = ',',
        header = TRUE,
        stringsAsFactors = FALSE
    )

    regionAreas <- NULL
    if( !is.null(regsFileName) & !is.null(gridFileName) )
        regionAreas <- computeRegionAreas(regsFileName, gridFileName)    

    param<-merge(omega, pntRate, by='region')
    
    if(!is.null(regionAreas))
        param<-merge(param, regionAreas, by='region')
    
    param<-merge(param, pop_reg, by='region')

    param[,dedupPntRate:=(omega1 + 0.5 * omega2) * pntRate]
    param[,alpha:=dedupPntRate * Nreg]
    param[,beta:=(1 - dedupPntRate) * Nreg]
    if(!is.null(regionAreas)) {
        param[,factor:= sqrt(1 + (2 * cv / (1 + rel_bias))^2) - 1]    
        param[,theta:= 0.5 * (Nreg / regionArea_km2) * (1 + rel_bias) * factor]
        param[,zeta := 2 / factor]    
        param[,factor:= NULL]
        param[,Q := regionArea_km2 * theta / (1 + regionArea_km2 * theta)]
    }
    return (param)
}