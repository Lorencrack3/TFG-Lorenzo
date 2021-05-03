#' @title Computes the origin-destination matrices.
#'
#' @description Computes the origin-destination matrices for all pairs of time instants time_from-time_to.
#'
#' @param nt0 The population at t0.
#'
#' @param nnetODFileName the name of the file where the population moving from one region to another is stored. This is
#'   an output of the \code{aggregation} package.
#'
#' @param zip If TRUE the file where where the population moving from one region to another is stored is a zipped csv
#'   file, otherwise it is  simple csv file.
#'
#' @param rndVal Controls if the random values generated for each t are returned or not in the result of this function.
#'   If TRUE, the random values generated according to the corresponding distribution are returned in the results, if
#'   FALSE only the summary statistics for each t and region are returned.
#'
#' @param ciprob Value of probability of the CI (between 0 and 1) to be estimated. If NULL the default value is 0.89.
#'
#' @param method The method to compute credible intervals. It could have 2 values, 'ETI' or 'HDI'. The default value is
#'   'ETI.
#'
#' @return A list with one element for each pair of time_from-time_to. Each element of the list is also a list with one
#'   or two elements, depending on the value of the rndVal parameter. If rndVal is TRUE there are two elements in the
#'   list corresponding to time instant a pair time_from-time_to. The first one is a data.table object with some
#'   descriptive statistics for the origin-destination matrix, containing the following columns:\code{region_from,
#'   region_to, Mean, Mode, Median, Min, Max, Q1, Q3, IQR, SD, CV, CI_LOW, CI_HIGH}. The second one is a data.table
#'   object with the random values for origin-destination matrix generated for each pair of time instants
#'   time_from-time_to and each pair of regions region_from-region_to, with the following columns: \code{region_from,
#'   region_to,iter, NPop}. If rndVal is FALSE the list for a pair of time instants time-from-time_to contains only the
#'   first element previously mentioned. The name of the list element corresponding to a pair of time instants is
#'   'time_from-time_to' and the name of the two list elements giving the descriptive statistics and random values are
#'   'stats' and 'rnd_values'.
#'
#' @import data.table
#' @import utils
#' @include computeStats.R
#' @include computeTau.R
#' @export
computePopulationOD <- function(nt0, nnetODFileName, zip = TRUE, rndVal = FALSE, ciprob = NULL, method = 'ETI') {
    
    tau = computeTau(nnetODFileName, zip)
    regs<-unique(nt0$region)
    n<-length(regs)
    
    nt0[,iter:=rep(1:(nrow(nt0)/n),times = n)]
    runningN <- copy(nt0)[, N := NULL]
    setcolorder(runningN, c('region', 'iter'))
    
    times_from <- sort(unique(tau$time_from))
    times_to <- sort(unique(tau$time_to))
    
    result <- vector(mode = 'list', length = length(times_from))
    
    names(result)<-paste0(times_from, "-",times_to)
    k=1
    pb <- txtProgressBar(min = 0, max = length(times_from), style = 3)
    for(timeFrom in times_from) {
        tempDT <- merge(tau[time_from == timeFrom], runningN,
                        by.x = c('region_from', 'iter'), by.y = c('region', 'iter'), 
                        allow.cartesian = TRUE)
        
        tempOD <- tempDT[, OD := tau_Nnet * NPop]

        runningN <- tempDT[, list(NPop = round(sum(tau_Nnet * NPop))), by = c('region_to', 'iter')]
        setnames(runningN, c('region_to'), c('region'))
        rm(tempDT)
        
        gc()
        tempOD<-tempOD[, NPop:=NULL]
        tempOD<-tempOD[, time_from:=NULL]
        tempOD<-tempOD[, time_to:=NULL]
        tempOD<-tempOD[, tau_Nnet:=NULL]
        setnames(tempOD, 'OD', 'NPop')
        stats_line <- computeStats(tempOD, ciprob, method,  by = c('region_from', 'region_to'))
        result[[k]][[1]] <- stats_line
        names(result[[k]]) <- 'stats'
        if(rndVal == TRUE) {
            setcolorder(tempOD, c('region_from', 'region_to', 'iter', 'NPop'))
            result[[k]][[2]] <- tempOD
            names(result[[k]]) <-c('stats', 'rnd_values')
        }
        setTxtProgressBar(pb, k)
        k <- k + 1
    }
    close(pb)
    return (result)
    
}