#' @title Computes the duplicity probabilities for each device using a Bayesian approach.
#'
#' @description Computes the duplicity probabilities for each device using a Bayesian approach. It uses two methods:
#'   "pairs" and "1to1". The "pairs" method considers the possible pairs of two compatible devices. These devices were
#'   selected by \code{computePairs()} function taking into consideration the antennas where the devices are connected
#'   and the coverage areas of antennas. Two devices are considered compatible if they are connected to the same or to
#'   neighbouring antennas. Thus, the data set with pairs of devices will be considerable smaller than all possible
#'   combinations of two devices. The "1to1" method considers all pairs of two devices when computing the duplicity
#'   probability, the time complexity being much greater than that of the "pairs" method. Both methods uses parallel
#'   computations to speed up the execution. They build a cluster of working nodes and splits the pairs of devices
#'   equally among these nodes.
#'
#' @param method Selects a method to compute the duplicity probabilities. It could have one of the two values: "pairs"
#'   or "1to1". When selecting "pairs" method, the pairs4dupl parameter contains only the compatible pairs of devices,
#'   i.e. the pairs that most of the time are connected to the same or to neighbouring antennas. "1to1" method checks
#'   all possible combinations between devices to compute the duplicity probabilities.
#'
#'
#' @param deviceIDs A vector with the device IDs. It is obtained by calling the \code{getDevices()} function.
#'
#' @param pairs4dupl A data.table object with pairs of devices and pairs of antennas where these devices are connected.
#'   It can be obtained by calling \code{computePairs()} function.
#'
#' @param modeljoin The joint HMM model returned by \code{getJointModel()} function.

#' @param llik A vector with the values of the log likelihoods after the individual HMM models for each device were
#'   fitted. This vector can be obtained by calling \code{fitModels()} function.
#'
#' @param P1 The apriori duplicity probility as it is returned by \code{aprioriDuplicityProb()} function. It is used
#'   when "pairs" method is selected.
#'
#' @param Pii Apriori probability of a device to be in a 1-to-1 correspondence with the holder as it is returned by
#'   \code{aprioriOneDeviceProb()} function. This parameter is used only when "1to1" method is selected.
#'
#' @param init A logical value. If TRUE, the \code{fit()} function uses the stored steady state as fixed initialization,
#'   otherwise the steady state is computed at every call of \code{fit()} function.
#'
#' @param lambda It is used only when "1to1" method is selected and a non NULL value mean that the computation of the
#'   duplicity probablities is performed according to the method described in \emph{An end-to-end statistical process
#'   with mobile network data for Official Statistics} paper.
#'
#' @return a data.table object with two columns: 'deviceID' and 'dupP'. On the first column there are deviceIDs and on
#'   the second column the corresponding duplicity probability, i.e. the probability that a device is in a 2-to-1
#'   correspondence with the holder.
#'
#' @import data.table
#' @import destim
#' @import parallel
#' @import Rsolnp
#' @export
computeDuplicityBayesian <-
  function(method,
           deviceIDs,
           pairs4dupl,
           modeljoin,
           llik,
           P1 = NULL,
           Pii = NULL,
           init = TRUE,
           lambda = NULL) {
    ndevices <- length(deviceIDs)
    jointEmissions <- emissions(modeljoin)
    noEvents <- apply(jointEmissions, 2, sum)
    colNamesEmissions <- colnames(jointEmissions)
    eee0 <- 1:length(colNamesEmissions)
    envEmissions <- new.env(hash = TRUE)
    
    for (i in eee0)
      envEmissions[[colNamesEmissions[i]]] <- i
    rm(colNamesEmissions)
    keepCols <- names(pairs4dupl)[-which(names(pairs4dupl) %in% c("index.x", "index.y"))]
    
    if (method == "pairs") {
      P2 <- 1 - P1
      alpha <- P2 / P1
      cl <-
        buildCluster(
          c(
            'pairs4dupl',
            'keepCols',
            'noEvents',
            'modeljoin',
            'envEmissions',
            'alpha',
            'llik',
            'init'
          ),
          env = environment()
        )
      ichunks <- clusterSplit(cl, 1:nrow(pairs4dupl))
      res <-
        clusterApplyLB(
          cl,
          ichunks,
          doPair,
          pairs4dupl,
          keepCols,
          noEvents,
          modeljoin,
          envEmissions,
          alpha,
          llik,
          init
        )
      stopCluster(cl)
      dupP.dt <- buildDuplicityTablePairs(res, deviceIDs)
    }
    
    else if (method == "1to1") {
      cl <-
        buildCluster(
          c(
            'pairs4dupl',
            'deviceIDs',
            'keepCols',
            'noEvents',
            'modeljoin',
            'envEmissions',
            'llik',
            'init'
          ),
          env = environment()
        )
      ichunks <- clusterSplit(cl, 1:ndevices)
      res <-
        clusterApplyLB(
          cl,
          ichunks,
          do1to1,
          pairs4dupl,
          deviceIDs,
          keepCols,
          noEvents,
          modeljoin,
          envEmissions,
          llik,
          init
        )
      stopCluster(cl)
      dupP.dt <- buildDuplicityTable1to1(res, deviceIDs, Pii, lambda)
    } else {
      stop("Method unknown!")
    }
    return(dupP.dt)
  }


do1to1 <-
  function(ichunks,
           pairs,
           devices,
           keepCols,
           noEvents,
           modeljoin,
           envEms,
           llik,
           init) {
    ndev <- length(devices)
    n <- length(ichunks)
    ll.matrix <- matrix(0L, nrow = n, ncol = ndev)
    for (i in 1:n) {
      ii <- ichunks[[i]]
      llik_ii <- llik[ii]
      if ((ii + 1) <= ndev)
        j_list <- (ii + 1):ndev
      else
        j_list <- c()
      for (j in j_list) {
        newevents <- sapply(pairs[index.x == ii & index.y == j, ..keepCols],
                            function(x)
                              ifelse(!is.na(x), envEms[[x]], NA))
        
        if (all(is.na(newevents)) |
            any(noEvents[newevents[!is.na(newevents)]] == 0)) {
          llij <- Inf
        } else {
          fitTry <-
            try(modeljoin_ij <-
                  fit(modeljoin,
                      newevents,
                      init = init,
                      method = "solnp"))
          if (inherits(fitTry, "try-error")) {
            llij <- Inf
          } else {
            llij <- logLik(modeljoin_ij, newevents)
          }
        }
        ll.aux_ij <- llik_ii + llik[j] - llij
        ll.matrix[i, j] <- ll.aux_ij
      } #end for j
    }
    return(ll.matrix)
  }

doPair <-
  function(ichunks,
           pairs,
           keepcols,
           noEvents,
           modeljoin,
           envEms,
           alpha,
           llik,
           init) {
    localdup <- copy(pairs[ichunks, 1:2])
    for (i in 1:length(ichunks)) {
      index.x0 <- localdup[i, index.x]
      index.y0 <- localdup[i, index.y]
      newevents <-
        sapply(pairs[ichunks[[i]], ..keepcols], function(x)
          ifelse(!is.na(x), envEms[[x]], NA))
      if (all(is.na(newevents)) |
          any(noEvents[newevents[!is.na(newevents)]] == 0)) {
        llij <- Inf
      } else {
        fitTry <-
          try(modeljoin_ij <-
                fit(modeljoin,
                    newevents,
                    init = init,
                    method = "solnp"))
        if (inherits(fitTry, "try-error")) {
          llij <- Inf
        } else {
          llij <- logLik(modeljoin_ij, newevents)
        }
      }
      dupP0 <-
        1 / (1 + (alpha * exp(llij - llik[index.x0] - llik[index.y0])))
      localdup[index.x == index.x0 &
                 index.y == index.y0, dupP := dupP0]
    }
    return(localdup)
  }
