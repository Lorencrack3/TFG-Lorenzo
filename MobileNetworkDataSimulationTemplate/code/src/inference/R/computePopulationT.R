#' @title Computes population counts at time  instants t >t0.
#'
#' @description Computes the distribution of the population counts for all times instants t > t0.
#'
#' @param nt0 The population at t0.
#'
#' @param nnetODFileName the name of the file where the population moving from one region to another is stored. This is
#'   an output of the \code{aggregation} package.
#'
#' @param zip If TRUE the file where where the population moving from one region to another is stored is a zipped csv
#'   file, otherwise it is  simple csv file.
#'
#' @param rndVal Controls if the random values generated for each t >t0 are returned or not in the result of this
#'   function. If TRUE, the random values generated according to the corresponding distribution are returned in the
#'   results, if FALSE only the summary statistics for each t>t0 and region are returned.
#'
#' @param ciprob Value of probability of the CI (between 0 and 1) to be estimated. If NULL the default value is 0.89.
#'
#' @param method The method to compute credible intervals. It could have 2 values, 'ETI' or 'HDI'. The default value is
#'   'ETI.
#'
#' @return A list with one element for each time instant (including t0). Each element of the list is also a list with
#'   one or two elements, depending on the value of the rndVal parameter. If rndVal is TRUE there are two elements in
#'   the list corresponding to time instant t. The first one is a data.table object with some descriptive statistics for
#'   the population count at time t, containing the following columns:\code{region, Mean, Mode, Median, Min, Max, Q1,
#'   Q3, IQR, SD, CV, CI_LOW, CI_HIGH}. The second one is a data.table object with the random values for population
#'   count generated for each region, with the following columns: \code{region, iter, NPop}. If rndVal is FALSE the list
#'   for time instant t contains only the first element previously mentioned. The name of the list element corresponding
#'   to time instant t is 't' and the name of the two list elements giving the descriptive statistics and random values
#'   for time t are 'stats' and 'rnd_values'.
#'
#' @import data.table
#' @import utils
#' @include computeStats.R
#' @include computeTau.R
#' @export
computePopulationT <- function(nt0, nnetODFileName, zip = TRUE, rndVal = FALSE, ciprob = NULL, method = 'ETI') {
    
    tau = computeTau(nnetODFileName, zip)
    regs<-unique(nt0$region)
    n<-length(regs)
    
    nt0[,iter:=rep(1:(nrow(nt0)/n),times = n)]
    
    runningN <- copy(nt0)[, N := NULL]
    setcolorder(runningN, c('region', 'iter'))
    
    times <- sort(unique(tau$time_from))
    timeIncrement <- unique(diff(times))
    
   
    result <- vector(mode = 'list', length = length(times)+1)
    names(result)<-as.character(c(times, times[length(times)]+timeIncrement))
    
    result[[1]][[1]] <- computeStats(runningN, ciprob, method)
    names(result[[1]]) <-c('stats')
    if(rndVal==TRUE) {
        result[[1]][[2]] <- runningN
        names(result[[1]]) <-c('stats', 'rnd_values')
    }


    k=1
    pb <- txtProgressBar(min = 0, max = length(times), style = 3)
    for(timeFrom in times[k:(length(times))]) {
        tempDT <- merge(tau[time_from == timeFrom], runningN,
                        by.x = c('region_from', 'iter'), by.y = c('region', 'iter'),
                        allow.cartesian = TRUE)
        runningN <- tempDT[, list(NPop = round(sum(tau_Nnet * NPop))), by = c('region_to', 'iter')]
        setnames(runningN, c('region_to'), c('region'))
        rm(tempDT)
        
        gc()
        
        stats_line <- computeStats(runningN, ciprob, method)
        k <- k + 1
        result[[k]][[1]] <- stats_line
        names(result[[k]]) <- 'stats'
        if(rndVal == TRUE) {
            result[[k]][[2]] <- runningN
            names(result[[k]]) <-c('stats', 'rnd_values')
        }
        setTxtProgressBar(pb, k)
    }
    close(pb)
    return (result)
}

doPopT <- function(ichunks, times, timeIncrement, tau, rndVal, ciprob, method ) {
    res<-vector(mode = 'list', length = length(ichunks))
    names(res) <- as.character(ichunks + timeIncrement)
    k <- 1            
    for(timeFrom in ichunks) {
        cat(paste0('Time ', timeFrom, '... '))  
        tempDT <- merge(tau[time_from == timeFrom], runningN,
                        by.x = c('region_from', 'iter'), by.y = c('region', 'iter'), 
                        allow.cartesian = TRUE)
        runningN <- tempDT[, list(NPop = sum(tau_Nnet * NPop)), by = c('region_to', 'iter')]
        setnames(runningN, c('region_to'), c('region'))
        rm(tempDT)
        gc()
        cat(paste0(' ok.\n'))
        
        stats_line <- computeStats(runningN, ciprob, method)
        res[[k]][[1]] <- stats_line
        names(res[[k]]) <- 'stats'

        if(rndVal == TRUE) {
            res[[k]][[2]] <- runningN
            names(res[[k]]) <-c('stats', 'rnd_values')
        }
        k <- k + 1
        
    }
    return(res)
}