#' @title Computes the distribution of the population count at initial time instant.
#'
#'
#' @description Computes the distribution of the population count at initial time instant using one of the three
#'   distributions: Negative Binomial, Beta Negative Binomial or State Process Negative Binomial.
#'
#'
#' @param nnet The random values generated with \code{aggregation} package for the number of individuals detected by the
#'   network.
#'
#' @param params The parameters of the distribution. It should be a data.table object with the following columns:
#'   \code{region, omega1, omega2, pnrRate, regionArea_km2, N0, dedupPntRate, alpha, beta, theta, zeta, Q}.
#'
#' @param popDistr The distribution to be used for population count. This parameter could have one of the following
#'   values: \code{NegBin} (negative binomial distribution), \code{BetaNegBin} (beta negative binomial distribution) or
#'   \code{STNegBin}  (state process negative binomial distribution).
#'
#' @param rndVal If FALSE the result return by this function will be a list with a single element, a data.table object
#'   with the following columns: \code{region, Mean, Mode, Median, SD, Min, Max, Q1, Q3, IQR, CV, CI_LOW, CI_HIGH}. If
#'   TRUE the list will have a second element which is a data.table object containing the random values generated for
#'   each region.
#'
#' @param ciprob Value of probability of the CI (between 0 and 1) to be estimated. If NULL the default value is 0.89.
#'
#' @param method The method to compute credible intervals. It could have 2 values, 'ETI' or 'HDI'. The default value is
#'   'ETI.
#'
#' @return A list object with one or two elements. If rndVal is FALSE the list will have a single element with
#'   descriptive statistics for the population count, which is a data.table object with the following columns:
#'   \code{region, Mean, Mode, Median, Min, Max, Q1, Q3, IQR, SD, CV, CI_LOW, CI_HIGH}. If rndVal is TRUE the list will
#'   have a second element which is a data.table object containing the random values generated for each region. The name
#'   of the two list elements giving the descriptive statistics and random values for time t are 'stats' and
#'   'rnd_values'.
#'
#' @import data.table
#' @import extraDistr
#' @import bayestestR
#' @include utils.R
#' @export
computeInitialPopulation <- function(nnet, params, popDistr, rndVal = FALSE, ciprob = NULL, method = 'ETI') {
    
    if(!( popDistr %in% c('NegBin', 'BetaNegBin', 'STNegBin')))
        stop('popDistr should have one of the following values: NegBin, BetaNegBin, STNegBin!')
        
    Ninit<-merge(nnet, params, by='region', all.x = TRUE, allow.cartesian = TRUE)
    
    if(popDistr == 'STNegBin') {
        NIP <- copy(Ninit)[, row := .I][
            , list(region = region,
                   N = N,
                   NPop = N + rnbinom(1, N + zeta + 1, 1 - beta * Q / (alpha + beta))), by = 'row']
    }
    
    if(popDistr == 'NegBin') {
        NIP <- copy(Ninit)[, row := .I][
            , list(region = region,
                   N = N,
                   NPop = N + rnbinom(1, N + 1, (alpha - 1) / (alpha + beta - 1))), by = 'row']
        
    }
    if(popDistr == 'BetaNegBin') {
        NIP <- copy(Ninit)[, row := .I][
            , list(region = region,
                   N = N,
                   NPop = N + rbnbinom(1, N + 1, alpha - 1, beta)), by = 'row']
        
    }
    NIP[,row:=NULL]    
    
    stats<-computeStats(NIP, ciprob, method)
    result<-list()
    result[[1]]<-stats
    names(result) <-'stats'
    if(rndVal) {
        result[[2]]<-NIP
        names(result) <-c('stats', 'rnd_values')
    }
    return(result)
    
}