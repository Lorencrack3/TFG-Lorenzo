#' @title  Builds a data.table object with device IDs and the duplicity probability for each device.
#'
#' @description Builds a data.table object with two columns: device IDs and the duplicity probability for each device.
#'   This function is called in case of using the "pairs" method of computing the duplicity probabilities, receiving a
#'   list of data.tables objects returned by each worker node with the device IDs on the first two columns and the
#'   corresponding duplicity probability on the third column. It is a utility function and it is not accesible from
#'   outside the package.
#'
#' @param res A list of data.table objects returned by each worker node with the device IDs on the first two columns and
#'   the corresponding duplicity probability on the third column for (the pair of devices in the first two columns).
#'
#' @param devices A vector with device IDs.
#'
#' @return A data.table object with two columns: the device IDs and the corresponding duplicity probability for each
#'   device.
buildDuplicityTablePairs <- function(res, devices) {
  dup <- NULL
  for (i in 1:length(res))
    dup <- rbind(dup, res[[i]])
  rm(res)
  
  dup[, deviceID1 := devices[index.x]][, deviceID2 := devices[index.y]]
  dup2 <- copy(dup)
  dup2[, deviceID1 := devices[index.y]][, deviceID2 := devices[index.x]]
  dup <- dup[, .(deviceID1, deviceID2, dupP)]
  dup2 <- dup2[, .(deviceID1, deviceID2, dupP)]
  dupProb.dt <- rbindlist(list(dup, dup2))
  
  allPairs <- expand.grid(devices, devices)
  rows_to_keep <- allPairs[, 1] != allPairs[, 2]
  allPairs <- as.data.table(allPairs[rows_to_keep,])
  setnames(allPairs, c(c("deviceID1", "deviceID2")))
  allDupProb.dt <-
    merge(allPairs[, .(deviceID1, deviceID2)],
          dupProb.dt,
          all.x = TRUE,
          by = c("deviceID1", "deviceID2"))
  rm(dup, dup2, dupProb.dt)
  allDupProb.dt[is.na(dupP), dupP := 0]
  dupProb.dt <- copy(allDupProb.dt)[, max(dupP), by = "deviceID1"]
  setnames(dupProb.dt, c("deviceID1", "V1"), c("deviceID", "dupP"))
  return(dupProb.dt)
  
}