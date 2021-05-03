#' @title Generates random values according to a Poisson multinomial probability
#'   distribution.
#'
#' @description Generates random values according to a Poisson multinomial
#'   probability distribution. A point estimation derived from this distribution
#'   (mean, mode) represents an estimation of the number of individuals detected
#'   by the network in a region. Regions are composed as a number of adjacent
#'   tiles. This is the only one function of this package available to users to
#'   compute an estimation of the number of detected individuals.
#'
#' @param n The number of random values to be generated.
#'
#' @param gridFileName The name of the .csv file with the grid parameters.
#'
#' @param dupFileName The name of the .csv file with the duplicity probability
#'   for each device. This is an output of the \code{deduplication} package.
#'
#' @param regsFileName The name of the .csv file defining the regions. It has
#'   two columns: \code{ tile, region}. The first column contains the IDs of
#'   each tile in the grid while the second contains the number of a region.
#'   This file is defined by the user and it can be created with any text
#'   editor.
#'
#' @param postLocPath The path where the files with the posterior location
#'   probabilities for each device can be found. A file with the location
#'   probabilities should have the name \code{prefix_ID.csv} where \code{ID} is
#'   replaced with the device ID and \code{prefix} is given as a parameter to
#'   this function.
#'
#' @param prefix A prefix that is used to compose the file name with posterior
#'   location probabilities.
#'
#' @param times A vector with the time instants when the events were
#'   registered. 
#'
#' @param seed The value of the random seed to be used by the random number generator.
#'
#' @return A data.table object with the following columns: \code{time, region,
#'   N, iter}. The last column contains the index of the random value (given in
#'   column \code{N}) generated for each time instant and region.
#'
#' @import data.table
#' @import deduplication
#' @import Matrix
#' @import  parallel
#' @include doAggr.R
#' @include buildCluster.R
#' @export
rNnetEvent <- function(n, gridFileName, dupFileName, regsFileName, postLocPath, prefix, times, seed = 123) {
  # 1. read duplicity probs.

  if (!file.exists(gridFileName))
    stop(paste0(gridFileName, " does not exists!"))
  
  gridParams <- readGridParams(gridFileName)
  
  if (!file.exists(dupFileName))
    stop(paste0(dupFileName, " does not exists!"))
  
  dupProbs <- fread(
    dupFileName,
    sep = ',',
    header = FALSE,
    stringsAsFactors = FALSE
  )
  setnames(dupProbs, c('device', 'dupP'))
  devices<-sort(as.numeric(dupProbs[,1][[1]]))
  
  # 2. read regions
  if (!file.exists(regsFileName))
    stop(paste0(regsFileName, " does not exists!"))
  
  regions <- fread(
    regsFileName,
    sep = ',',
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # 3. read posterior location probabilities
  ndevices <- nrow(dupProbs)
  postLoc <- NULL
  ntiles <- gridParams$nrow * gridParams$ncol
  ntimes <- length(times)
  timeincr<-times[2]-times[1]
  prefix <- prefix
  postLocPath <-postLocPath
  
  cl <- buildCluster( c('postLocPath', 'prefix', 'devices', 'ntiles', 'ntimes', 'timeincr') , env=environment() )
  ich<-clusterSplit(cl, 1:ndevices)

  res<-clusterApplyLB(cl, ich, doRead1, postLocPath, prefix, devices, ntiles, ntimes, timeincr)

  for(i in 1:length(res)) {
    postLoc <- c(postLoc, unlist(res[[i]]))
  }

  T <- ncol(postLoc[[1]])
   if(T != length(times))
     stop("Inconsistent data provided: the length of times vector is not the same as the number of time instants computed from the posterior location probabilities files")

  tiles <- 0:(ntiles-1)
  # 4. computation begins ...
  clusterExport(cl,  c('dupProbs', 'regions', 'postLoc', 'tiles'), envir=environment())
  clusterSetRNGStream(cl, iseed=seed)
  ichunks <- clusterSplit(cl, 1:T)
#return(list(cl, ichunks, doAggr, n, tiles, postLoc, dupProbs, regions))
  res <-
     clusterApplyLB(
       cl,
       ichunks,
       doAggr,
       n,
       tiles,
       postLoc,
       dupProbs,
       regions
     )
  stopCluster(cl)

  result <- rbindlist(res)
  indices <- result$time
  result$time <- times[indices]
  return (result)
}


 doRead1 <- function(ichunks, postLocPath, prefix, devices, ntiles, ntimes, timeincr) {
    n <- length(ichunks)
    local_postLoc<-list(length = n)
    k<-1
   for( i in ichunks) {
      tmp <- readPostLocProb(postLocPath, prefix, devices[i])
      local_postLoc[[k]] <- sparseMatrix(i=tmp$tile, j=(tmp$time)/timeincr, x=tmp$probL, dims=c(ntiles, ntimes), index1=FALSE)
      k <- k + 1
    }
   return (local_postLoc)
 }
