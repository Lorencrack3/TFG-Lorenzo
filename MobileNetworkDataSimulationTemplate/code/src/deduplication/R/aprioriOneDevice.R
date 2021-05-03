#'@title Apriori probability of a device to be in a 1-to-1 correspondence with the holder.
#'
#'@description Apriori probability of a device to be in a 1-to-1 correspondence with the holder. It is computed simply
#'  as \eqn{(2*Ndev1-Ndev)/Ndev1} where \eqn{Ndev1} is the number of devices that are in a 1-to-1 correspondence with the persons
#'  that hold them and \eqn{Ndev} is the total number of devices.
#'
#'@param prob2Devices The probability of a person to have 2 devices. In case of using simulated data, this parameter is
#'  read from the simulation configuration file.
#'
#'@param ndevices The number of devices detected by the network during the time horizon under consideration.
#'
#'@return The apriori probability of a device to be in a 1-to-1 correspondence with the holder.
#'
#'@export
aprioriOneDeviceProb <- function(prob2Devices, ndevices) {
  n_ext <- round((1 - prob2Devices) * ndevices)
  Pii <- (2 * n_ext - ndevices) / (n_ext)  # apriori probability of 1:1
  return(Pii)
}