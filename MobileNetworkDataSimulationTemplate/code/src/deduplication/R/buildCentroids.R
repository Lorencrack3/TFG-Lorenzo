#' @title Builds the centroid of the rectangular tiles.
#'
#' @description Builds the centroid of each (rectangular) tile in the grid. The centroid is simple the center of the
#'   tile, i.e. the intersection poit of the diagonals of the tile.
#'
#' @param ntiles_x The number of tiles on OX.
#'
#' @param ntiles_y The number of tiles on OY.
#'
#' @return A data.table object with three columns: tile ID, centroidCoord_x and centroidCoord_y.
#'
#' @import data.table
#' @export
buildCentroids <-
  function(ntiles_x,
           ntiles_y,
           tile_sizeX,
           tile_sizeY) {
    
    ntiles <- ntiles_x * ntiles_y
    centroidCoord_x  <- (0.5 + 0:(ntiles_x - 1)) * tile_sizeX
    centroidCoord_y  <- (0.5 + 0:(ntiles_y - 1)) * tile_sizeY
    centroidCoord <-as.data.table(expand.grid(x = centroidCoord_x, y = centroidCoord_y))[, tile := 0:(ntiles - 1)]
    return(centroidCoord)
  }
