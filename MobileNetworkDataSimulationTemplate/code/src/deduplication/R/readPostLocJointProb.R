#' @title Read a file with the posterior location joint probabilities.
#'
#' @description Read a csv file with the posterior location joint probabilities.
#'  Each row of the file corresponds to a pair of tiles and a time instant.
#'
#' @param path The path to the location where the posterior location probabilities are stored. The file with the
#'   location probabilities should have the name \code{postLocDevice_ID.csv} where \code{ID} is replaced with the device
#'   ID.
#' @param prefixName The file name prefix. The whole file name is composed by a concatenation of \code{prefixName},
#'   \code{_} and \code{deviceID}.
#'   
#' @param deviceID The device ID for which the posterior location probabilities are read.
#'
#' @return A Matrix object with the posterior location probabilities for the device with ID equals to deviceID. A row
#'   corresponds to a tile and a column corresponds to a time instant.
#'
#' @import data.table
#' @export
readPostLocJointProb <-function(path, prefixName, deviceID) {
  file <- paste0(path,"/", prefixName, "_", as.character(deviceID), ".csv")
  if(!file.exists(file))
    stop (paste0('file with posterior location probabilities files does not exist ', file))
  
  postLoc <- fread(file, sep = ',',stringsAsFactors = FALSE,header = FALSE)
  setnames(postLoc, paste0('V', c(1:5, ncol(postLoc))), 
           c('device', 'time_from', 'time_to', 'tile_from', 'tile_to', 'probL'))
  
  return (postLoc)
  
}