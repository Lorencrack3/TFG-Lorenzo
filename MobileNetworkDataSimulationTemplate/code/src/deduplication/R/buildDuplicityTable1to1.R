#' @title  Builds a data.table object with device IDs and the duplicity probability for each device.
#'
#' @description Builds a data.table object with with two columns: device IDs and the duplicity probability for each
#'   device. This function is called in case of using the "1to1" method of computing the duplicity probabilities,
#'   receiving a list of matrix objects returned by each worker node with lines corresponding to a subset of devices.
#'   These parts a put together and form a symmetic matrix used to compute the duplicity probability. It is a utility
#'   function and it is not accesible from outside the package.
#'
#' @param res A list of matrix objects returned by each worker node with the number of rows corresponding to the number
#'   of devices allocated to an working node and the number of columns equals to the number of devices.
#'
#' @param devices A vector with device IDs.
#' 
#' @param Pii Apriori probability of devices to be in 1-to-1 correspondence with the owner.
#' 
#' @param lambda If this parameter is non NULL, the duplicity probabilities are computed according to the approach described in
#' "An end-to-end statistical process with mobile network data for Official Statistics" paper.
#'
#' @return A data.table object with two columns: the device IDs and the corresponding duplicity probability for each
#'   device.
#'   
buildDuplicityTable1to1 <- function(res, devices, Pii, lambda = NULL) {
  ndevices <- length(devices)
  
  if(!is.null(Pii)) {
    Pij <- (1 - Pii) / (ndevices - 1)    # priori prob. of duplicity 2:1
    alpha <- Pij / Pii
  } else {
    if(!is.null(lambda)) {
      alpha <- vector(length = ndevices)
      if(length(lambda) == 1) {
        alpha <- rep(1/(lambda *(ndevices-1)), times = ndevices)
      }
      else {
        alpha <- 1/(lambda *(ndevices-1))
      }
    }
  }
  
  matsim <- NULL
  for (i in 1:length(res)) {
    matsim <- rbind(matsim, res[[i]])
  }
  rm(res)
  
  matsim[lower.tri(matsim)] <- t(matsim)[lower.tri(matsim)]
  
  dupP.dt <- data.table(deviceID = devices, dupP = rep(0, ndevices))
  
  for (i in 1:ndevices) {
    ll.aux <- matsim[i, -i]
    if(!is.null(Pii)) {
      dupP.dt[deviceID == devices[i], dupP := 1 - 1 / (1 + (alpha * sum(exp(ll.aux))))]
    }
    else {
      alpha2<-alpha[-i]
      dupP.dt[deviceID == devices[i], dupP := 1 - 1 / (1 + (sum(alpha2 * exp(ll.aux))))]
    }
  }
  return (dupP.dt)
  
}