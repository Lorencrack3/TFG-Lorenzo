#' @title Generates a sequence of Poisson multinomial random values for a
#'   sequence of time instants.
#'
#' @description Generates a sequence of Poisson multinomial random values for a
#'   sequence of time instants. This function is called internally by
#'   \code{rNnetEvent} function which partitions the time instants into equal
#'   chunks, creates a cluster of working processes and distributes the chunks
#'   to different workers that run in parallel. This is an internal function
#'   that cannot be directly called by the users.
#'
#' @param ichunks A partition of time instants.
#'
#' @param  n The number of random values to be generated.
#'
#'
#' @param tiles The list of tile indexes.
#'
#' @param postLoc A list with posterior location probabilities matrices. Each
#'   element of the list corresponds to a device.
#'
#' @param dupProbs The duplicity probabilities for each device. It is a
#'   data.table object with two columns: \code{deviceID, dupProb}.
#'
#' @param regions a data.table object with two columns: \code{tile, region}. It
#'   defines each region as a composition of multiple tiles.
#'
#'
#' @return a data.table object with the following column: \code{time, region,
#'   N}. It contains \code{n} random values for each combination
#'   \code{time-region}. An estimation of the number of individuals in a region
#'   can be obtained using the mean, mode or median of these values.
#'
#'
#'
#' @import data.table
#' @import deduplication
#' @include rNnet_Event.R
doAggr <- function(ichunks, n, tiles, postLoc, dupProbs, regions) {
  nIndividualsT <- list(length=length(ichunks))
  k<-1
  nTiles <- length(tiles)
  devices<-sort(as.numeric(dupProbs[,1][[1]]))
  ndevices <- length(devices)
  for(t in ichunks) {
    dedupLoc <- data.table()
    for(j in 1:ndevices) {
      x <- cbind(tiles, postLoc[[j]][,t], rep(devices[j], times = nTiles))
      dedupLoc <- rbind(dedupLoc, x)
      x <- NULL
    }
    colnames(dedupLoc)<-c('tile', 'postL', 'device')
    dedupLoc2_1 <- merge(dedupLoc, dupProbs, by = 'device', all.x = TRUE)
    dedupLoc1_1 <- copy(dedupLoc2_1)[ , singleP := 1- dupP][, dupP := NULL]
    dedupLoc2_1[,prob := postL * dupP][, devCount := 0.5][, postL := NULL][,dupP := NULL]
    dedupLoc1_1[,prob := postL * singleP][, devCount := 1][, postL := NULL] [, singleP := NULL]
    dedupProbs <- rbind(dedupLoc1_1, dedupLoc2_1)
    rm(dedupLoc2_1, dedupLoc1_1)
    dedupProbs <- merge( dedupProbs, regions, by = c('tile'))
    dedupProbs <- dedupProbs[ ,list(prob = sum(prob)), by = c('device', 'region', 'devCount')]
    nIndividuals_MNO <- as.data.table(rNnet_Event(n, dedupProbs))
    nIndividuals_MNO_molten <- melt( nIndividuals_MNO, variable.name = 'region', value.name = 'N')
    nIndividuals_MNO_molten <- cbind(rep(t, times = nrow(nIndividuals_MNO_molten)), nIndividuals_MNO_molten)
    nIndividuals_MNO_molten <- cbind(nIndividuals_MNO_molten, 1:n)
    colnames(nIndividuals_MNO_molten)<-c('time', 'region', 'N', 'iter')
    nIndividualsT[[k]]  <- nIndividuals_MNO_molten
    k<-k+1
  }
  return(rbindlist(nIndividualsT) )
}
