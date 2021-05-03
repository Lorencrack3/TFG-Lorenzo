#' @title Apriori probability of duplicity.
#'
#' @description Apriori probability of duplicity, i.e. the probability of a device to be in a 2-to-1 correspondence with
#'   the person who holds it. It is computed as \eqn{ 2 * (Ndev-Ndev2) / (Ndev * (Ndev-1)) } where \eqn{Ndev} is the
#'   total number of devices and \eqn{Ndev2} is the number of devices that are in a 1-to-1 correspondence with the
#'   persons that hold them.
#'
#' @param prob2Devices The probability of a person to have 2 devices. In case of using simulated data, this parameter is
#'   read from the simulation configuration file.
#'
#' @param ndevices The number of devices detected by the network during the time horizon under consideration.
#'
#' @return Apriori probability of duplicity, i.e. the probability of a device to be in a 2-to-1 correspondence with the
#'   person who holds it.
#'
#' @export
aprioriDuplicityProb <- function(prob2Devices, ndevices) {
  n_ext <- round((1 - prob2Devices) * ndevices)
  P1 <- (2 * (ndevices - n_ext)) / (ndevices * (ndevices - 1))
  return (P1)
}
