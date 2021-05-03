#'@title Computes the duplicity probabilities for each device using the trajectory approach.
#'
#'@description  Computes the duplicity probabilities for each device using the trajectory approach described in
#'  \href{https://webgate.ec.europa.eu/fpfis/mwikis/essnetbigdata/images/f/fb/WPI_Deliverable_I3_A_proposed_production_framework_with_mobile_network_data_2020_05_31_draft.pdf}{WPI
#'   Deliverable 3}.
#'
#'@param path The path where the files with the posterior location probabilities for each device are to be found.
#'
#'@param devices A vector with device IDs.
#'
#'@param gridParams A list with the number of rows and columns of the grid and the tile dimensions on OX and OY axes.
#'  The items of the list are named 'nrow', 'ncol', 'tileX' and 'tileY'.
#'
#'@param pairs A data.table object containing pairs with the IDs of compatible devices. It is obtained calling
#'  \code{buildPairs()} function.
#'
#'@param P1 The apriori probabilitity for a device to be in 1-to-1 correspondence with its owner.
#'
#'@param T The sequence of time instants in the data set.
#'
#'@param gamma A coefficient needed to compute the duplicity probability. See
#'  \href{https://webgate.ec.europa.eu/fpfis/mwikis/essnetbigdata/images/f/fb/WPI_Deliverable_I3_A_proposed_production_framework_with_mobile_network_data_2020_05_31_draft.pdf}{WPI
#'   Deliverable 3}.
#'
#'@return  a data.table object with two columns: 'deviceID' and 'dupP'. On the first column there are deviceIDs and on
#'  the second column the corresponding duplicity probability, i.e. the probability that a device is in a 2-to-1
#'  correspondence with the holder.
#'
#'@include centerOfProbabilities.R
#'@include readPostLocProb.R
#'@include dispersionRadius.R
#'@include buildCentroidProbs.R
#'@include buildCentroids.R
#'@include buildDeltaProb.R
#'@import data.table
#'@import parallel
#'@export
computeDuplicityTrajectory <-function(path, prefix, devices, gridParams, pairs, P1 , T, gamma) {
  devices <- sort(as.numeric(devices))
  centrs <- buildCentroids(gridParams$ncol, gridParams$nrow, gridParams$tileX, gridParams$tileY)
  ndevices <- length(devices)
  
  postLoc <- NULL
  centerOfProbs <- NULL
  dr<-NULL
  cpp<-NULL
  
  P2 <- 1 - P1
  alpha<-P1/P2

  ### from this point the code should go parallel  
  path<-path
  prefix<-prefix
  timeincr<-T[[2]]-T[[1]]
  ntimes<-length(T)
  ntiles<-gridParams$nrow * gridParams$ncol
  cl <- buildCluster( c('path', 'prefix', 'devices', 'centrs', 'timeincr', 'ntimes', 'ntiles'), env = environment())
  ichunks <- clusterSplit(cl, 1:ndevices)
  res <- clusterApplyLB( cl, ichunks, doLocations, path, prefix, devices, centrs, timeincr, ntimes, ntiles )
  for(i in 1:length(res)) {
    postLoc <- c(postLoc, res[[i]]$postLoc)
    centerOfProbs <- c(centerOfProbs, res[[i]]$centerOfProbs)
    dr <- c(dr, res[[i]]$dr)
  }
  rm(res)

  ichunks2 <- clusterSplit(cl, 1:ntimes)
  clusterExport(cl, varlist = c('postLoc'), envir = environment())
  res <- clusterApplyLB( cl, ichunks2, doCPP, centrs, postLoc )
  for(i in 1:length(res)) {
    cpp <- c(cpp, res[[i]])
  }
  rm(res)


  pairs<-pairs[,1:2]
  ichunks3 <- clusterSplit(cl, 1:nrow(pairs))
  
  clusterExport(cl, varlist = c('pairs', 'cpp', 'dr', 'ndevices', 'gamma', 'alpha'), envir = environment())
  res<-clusterApplyLB(cl, ichunks3, doPairs, pairs, cpp, ndevices, dr, ntimes, alpha, gamma)
  stopCluster(cl)
  dup<-buildDuplicityTablePairs(res, devices)
  return (dup)
  
}

doLocations <- function(ichunks, path, prefix, devices, centrs, timeincr, ntimes, ntiles ) {
  n <- length(ichunks)
  local_postLoc<-list(length = n)
  local_centerOfProbs<-list(length = n)
  local_dr <- list(length = n)
  k<-1
  for( i in ichunks) {
    tmp <- readPostLocProb(path, prefix, devices[i])
    local_postLoc[[k]] <-  sparseMatrix(i=tmp$tile, j=(tmp$time)/timeincr, x=tmp$probL, dims=c(ntiles, ntimes), index1=FALSE)
    local_centerOfProbs[[k]] <- centerOfProbabilities(centrs, local_postLoc[[k]])
    local_dr[[k]] <- dispersionRadius(centrs, local_postLoc[[k]], local_centerOfProbs[[k]])
    k <- k + 1
  }
  res<-list(postLoc = local_postLoc, centerOfProbs = local_centerOfProbs, dr = local_dr )
  return (res)
}

doCPP <- function(ichunks, centrs, postLoc) {
  n <- length(ichunks)
  local_cpp <- list(length = n)
  k<-1
  for ( i in ichunks) {
     local_cpp[[k]]<-buildCentroidProbs(centrs, postLoc, i)
     k<-k+1
  }
  return (local_cpp)
}


doPairs <- function(ichunks, pairs, cpp, ndev, dr, ntimes, alpha, gamma) {
  n <- length(ichunks)  
  localdup <- copy(pairs[ichunks, 1:2])
  for( i in ichunks ) {
    index_i <-pairs[i,1][[1]]
    index_j <-pairs[i,2][[1]]
    s1<-0
    for(t in 1:ntimes) {
      mm<-gamma*max(dr[[index_i]][t], dr[[index_j]][t])
      mdelta<-buildDeltaProb(cpp[[t]][[index_i]], cpp[[t]][[index_j]])
      s1<-s1+(abs(modeDelta(mdelta[[1]]))<mm & abs(modeDelta(mdelta[[2]]))<mm)
    }
    tmp <- s1/ntimes
    dupP0 <- 1 - 1/(1+alpha*tmp/(1-tmp))
    localdup[index.x == index_i &index.y == index_j, dupP := dupP0]
  }
  return (localdup)
}
