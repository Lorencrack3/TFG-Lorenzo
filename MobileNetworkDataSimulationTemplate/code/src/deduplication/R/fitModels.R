#' @title Fits the HMM model for each device.
#'
#' @description Fits the HMM model for each device using the \code{fit()} function from \pkg{destim} package. The
#'   computations are done in parallel to reduce the running time using all the available cores. This function creates a
#'   cluster of working nodes, splits the devices equally among the working nodes and assigns a partition of devices to
#'   each node in the cluster. For Unix-like operating systems, this functions uses a "FORK" cluster while for Windows
#'   it uses a "SOCK" cluster.
#'
#' @param ndevices The number of devices.
#'
#' @param model The HMM model returned by \code{getGenericModel()} function. This model is fittted for each device.
#'
#' @param connections A matrix whose elements are the antenna ID where a device is connected at every time instant. This
#'   matrix is returned by \code{getConnections()} function.
#'
#'
#' @return A vector of log likelihoods computed using the fitted model for each device.

#' @import parallel
#' @import destim
#' @import Rsolnp
#' @export
fitModels <-
  function(ndevices, model, connections) {
      cl<-buildCluster(c('connections', 'model'), env = environment())
      ichunks <- clusterSplit(cl, 1:ndevices)
      res <- clusterApplyLB(cl, ichunks, doFit, model, connections)
      stopCluster(cl)
      return (unlist(res))
  }


doFit <- function(index, model, connections) {
  local_ll <- vector(length = length(index))
  k = 1
  for (j in index) {
    fitTry <-try(modeli <- fit(model,connections[j,], init = TRUE, method = "solnp"))
    if (inherits(fitTry, "try-error")) {
      local_ll[k] <- Inf
    } else {
      local_ll[k] <- logLik(modeli, connections[j, ])
    }
    k <- k + 1
  }
  return (local_ll)
}
