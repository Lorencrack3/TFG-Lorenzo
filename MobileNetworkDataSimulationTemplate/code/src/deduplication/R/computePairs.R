#' @title Builds pairs of devices and corresponding connections.
#'
#' @description Builds a data.table object that contains pairs of antenna IDs in the form "antennaID1-antennaID2", where
#'   antennaID1 corresponds to a device while antennaID2 corresponds to another device, for each time instant over a
#'   time period when the network events are registered. These pairs are build for each distinct combination
#'   "deviceID1-deviceID2" of devices. The rows of the data.table object corresponds to a combination of devices and the
#'   columns correponds to different time instants. The first two columns contains the device IDs of each pair of
#'   devices and the rest of the columns correspond to a time instant and contains the pairs antennaID1-antennaID2 where
#'   the two devices are connected at that time instant.
#'
#' @param connections A matrix with the antenna IDs where the mobile devices are connected at every time instant. Each
#'   row corresponds to a device and each column to a time instant. This matrix is obtained by calling
#'   \code{getConnections()} function.
#'
#' @param ndevices The number of devices registered by the network. A vector with device IDs can be obtained by calling
#'   \code{getDevices()} function and the number of devices is simply the lenght of this vector.
#'
#' @param oneToOne If TRUE, the result is built to apply the method "1to1" to compute the duplicity probability for each
#'   device. This means that the result will contain all combinations of devices. If FALSE, the result will consider the
#'   proximity of antennas through the parameter antennaNeighbors and remove all the combinations of devices that are
#'   impossible to belong to a single person. If most of the time instants two devices are connected to two antennas
#'   that are not neighbours (i.e. their cells don't overlap) we consider that these two devices belong to different
#'   persons and remove this combination of devices from the result. The term "most of the time instants" is implemented
#'   like this: we add how many times during the time horizon two devices are connected to neighbouring antennas and
#'   than keep in the final result only those combinations of devices with this number greater than the quantile of the
#'   sequence 0...(Number of time instants) with probability 1-P1-limit. In this way we reduce the time complexity of
#'   the duplicity probability computation.
#'
#' @param P1 The apriori probability of duplicity. It is obtained by calling \code{aprioriDuplicityProb()} function.
#'
#' @param limit A number that stands for the error in computing apriori probability of duplicity.
#'
#' @param antennaNeighbors A data.table object with a single column 'nei' that contains pairs of antenna IDs  that are
#'   neighbours in the form antennaID1-antennaID2. We consider that two antennas are neighbours if their coverage
#'   areas has a non void intersection.
#'
#' @return a data.table object. The first two columns contain the devices indices while the rest of the columns contains
#'   pairs antennaID1-antennaID2 with antenna IDs where the devices are connected. There is one column for each time
#'   instant for the whole time horizon.
#'
#'
#' @import data.table
#' @export
computePairs <- function(connections, ndevices, oneToOne = TRUE, P1 = 0, limit = 0.05, antennaNeighbors = NULL) {
  
  if(oneToOne == FALSE & is.null(antennaNeighbors))
    stop("To apply \"pairs\" method you need to provide antennaNeighbors parameter")
  
  pairs4duplicity <- NULL
  
  connections[is.na(connections)] <- "NA"
  connections <- data.table(connections)[, index := .I]
  
  # make combinations of all devices
  allPairs <- data.table(t(combn(c(1:ndevices), 2)))
  setnames(allPairs, c("index.x", "index.y"))
  
  #antennaId where the first device is connected
  allPairs_connections.dt1 <- merge(allPairs, connections , by.x = c("index.x"),
                                    by.y = c("index"), all.x = TRUE)
  
  #antennaId where the second device is connected
  allPairs_connections.dt2 <- merge(allPairs, connections , by.x = c("index.y"),
                                    by.y = c("index"), all.x = TRUE)
  
  #change the column names to allow rbindlist
  setcolorder(allPairs_connections.dt2, names(allPairs_connections.dt1))
  allPairs_connections.dt <- rbindlist(list(allPairs_connections.dt1, allPairs_connections.dt2))
  
  #no longer needed
  rm(allPairs_connections.dt1)
  rm(allPairs_connections.dt2)
  
  pairs4duplicity <- allPairs_connections.dt[, lapply(.SD, paste0, collapse="-"), by = c("index.x", "index.y")]
  pairs4duplicity[pairs4duplicity == "NA-NA"] <- NA

  if(oneToOne == FALSE) {
    keepCols <- names(pairs4duplicity)[-which(names(pairs4duplicity) %in% c("index.x", "index.y"))]
    number <- sapply(pairs4duplicity[, ..keepCols],function(x){ x %in% antennaNeighbors[,nei]})
    pairs4duplicity[, number := apply(number, 1, sum)]
    rm(number)
    
    pairs4duplicity[, num_NA := Reduce(`+`, lapply(.SD,function(x) is.na(x))), .SDcols = keepCols]
    pairs4duplicity[, final_num := number + num_NA]

    pairs4duplicity <- pairs4duplicity[final_num >= quantile(c(0:length(keepCols)), probs = 1-P1-limit)][, c( "index.x", "index.y", keepCols), with = FALSE]
  } 
  return (pairs4duplicity)
}
