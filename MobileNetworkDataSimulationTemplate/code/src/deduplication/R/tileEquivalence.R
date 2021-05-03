#' @title Transforms the tiles indices from the notation used by the simulation software to the one used by the raster
#'   package.
#'
#' @description In order to perform the population estimations, the area of interest is overlapped with a rectangular
#'   grid of tiles. Each tile is a rectangle with predefined dimensions. This function is a utility function which
#'   transform the tiles indices from the numbering system used by the simulation software to the one used by the
#'   \pkg{raster} package. The simulation software uses a notation where the tile with index 0 is the bottom left tile while
#'   the raster package uses another way to number the tiles, tiles being numbered starting with 1 for the upper left
#'   tile.
#'
#' @param nrow Number of rows in the grid overpalling the area of interest.
#'
#' @param ncol Number of columns in the grid overpalling the area of interest.
#'
#' @return Returns a data.frame object with two columns: on the first column are the tile indices according to the
#'   raster package numbering and on the second column are the equivalent tile indices according to the simulation
#'   software numbering.
#'
#'
#' @export
tileEquivalence <- function(nrows, ncols){
  
  if(nrows <=0 )
    stop("nrows should be a positive number")
  if(ncols <=0 )
    stop("ncols should be a positive number")
  
  nCells <- nrows * ncols
  
  tileID_raster <- as.vector(matrix(1:nCells, ncol = ncols, byrow = TRUE))
  tileID_simulator <- as.vector(matrix(0:(nCells - 1), ncol = ncols, byrow = TRUE)[nrows:1, ])
  
  tileCorresp <- cbind(rasterCellID = tileID_raster, tileID = tileID_simulator)
  tileCorresp <- tileCorresp[order(tileCorresp[, 'rasterCellID']), ]
  colnames(tileCorresp)<-c('rasterCell', 'tile')
  return(tileCorresp)
}
