#' @title Calculate the distance between two cover areas in terms of tiles.
#'
#' @description This function obtains the number of tiles that are between two
#' cover areas by taking into account that a tile has signal from an antenna
#' in the center of the tile has signal.
#' 
#' @param from the antenna from which the device is going. It has to correspond to 
#' the name of a column in emissionProbs.
#'
#' @param to the antenna to which the device is going. It has to correspond to 
#' the name of a column in emissionProbs.
#' 
#' @param emissionProbs the matrix of emission probabilities.
#' 
#' @param centroidCoord.dt a data.table with the information about the coordinates
#' of the center of the tiles (rastercells). The required columns are: rasterCell,
#' centroidCoord_x and centroidCoord_y.
#' 
#' @param tile_size the size of each tile in meters.
#'
#' @return A numeric value with the number of tiles between the two cover areas.
#'
#' @examples
#' \dontrun{
#' }
#'
#' @import data.table 
#'
#'
#' @export
distance_coverArea <- function(from, to, emissionProbs, centroidCoord.dt, tile_size){


  # Check that neither from nor to are NA
  if(is.na(from) | is.na(to)){
    
    stop(paste0('[distance_coverArea] The arguments from or to cannot be NA.\n'))

  }
  # Check if from and to are in emissionProbs
  columnsLack <- setdiff(c(from, to), colnames(emissionProbs))
  if (length(columnsLack) > 0) {
    
    stop(paste0('[distance_coverArea] The following observedValues are missing in emissionProbs: ', 
                paste0(columnsLack, collapse = ', '),  '.\n'))
    
  }
  # Check if emissionProbs meets the restrictions
  if (!all(round(apply(emissionProbs, 1, sum), 6)==1 || round(apply(emissionProbs, 1, sum), 6)==0)) {
    
    stop(paste0('[distance_coverArea] There is some problem of consistency in emissionProbs\n'))
  }
  
  # Check if centroidCoord.dt has the required columns
  if(!all(c("rasterCell", "centroidCoord_x", "centroidCoord_y") %in% colnames(centroidCoord.dt))){
    
    stop(paste0('[distance_coverArea] centroidCoord.dt must have the columns: rasterCell, centroidCoord_x, centroidCoord_y.\n'))
    
  }  
  
  rc_from <- which(emissionProbs[, as.character(from)] > 0)
  rc_to <- which(emissionProbs[, as.character(to)] > 0)
  
  
  coord_from <- centroidCoord.dt[rasterCell %in% rc_from, 
                                 .(centroidCoord_x, centroidCoord_y)]
  coord_to <- centroidCoord.dt[rasterCell %in% rc_to, 
                                 .(centroidCoord_x, centroidCoord_y)]
  coords.list <- list(coord_from, coord_to)
  shorter.index <- which.min(c(nrow(coord_from), nrow(coord_to)))
  larger.index <- ifelse(shorter.index == 1, 2, 1)
  dista <- numeric(length = nrow(coords.list[[shorter.index]]))
  for(i in 1:nrow(coords.list[[shorter.index]])){
    
    D <- as.matrix(dist(rbind(coords.list[[shorter.index]][i,], 
                              coords.list[[larger.index]])))
    dista[i] <- min(D[-1, 1])
    
  }
  
  dist_tiles <- floor((min(dista) - tile_size)/tile_size)
  
  return(dist_tiles)

}