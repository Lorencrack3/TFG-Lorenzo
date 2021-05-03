#' @title Builds a matrix object containing the IDs of the antennas to which devices are connected.
#'
#' @description Builds a matrix object containing the IDs of the antennas to which devices are connected.
#' The number of rows equals the number of devices and the number of columns equals the number of time instants
#' when the network events were recorded. An element \code{[i,j]} in the returned matrix equals the ID of the antenna
#' where the mobile device with index \code{i} in the ordered list of device IDs (returned by \code{getDeviceIDs()}) is
#' connected at the time instant with index \code{j} in the sequence of time instants when the network events
#' were recorded.
#'
#' @param events A data.table object returned by \code{readEvents()} function.
#'
#' @return A matrix object with the antenna IDs where devices are connected for every time instants in the events file.
#' If a device is not connected to any antenna at a time instant, the corresponding element in the matrix will have the
#' value NA.
#' 
#'
#' @import data.table
#' @include getDeviceIDs.R
#' @export
getConnections <- function(events) {
  times <- unique(events[, 'time'])
  times <- sort(unlist(times[,'time']))
  
  devices <- getDeviceIDs(events)

  n1 <- length(times)
  n2 <- length(devices)
  
  connections <- matrix(-1L, ncol = n1, nrow = n2)
  
  evv <- matrix(ncol = 2, nrow = n1)
  evv[, 1] <- times
  colnames(evv) <- c('time', 'antennaID')
  
  for (i in 1:n2) {
    evv <-
      merge(evv, events[deviceID == devices[i], ][, c(1, 4)], by = 'time', all.x = TRUE)[, c(1, 3)]
    connections[i, ] <- t(evv[, 2])
  }

    return (connections)
}