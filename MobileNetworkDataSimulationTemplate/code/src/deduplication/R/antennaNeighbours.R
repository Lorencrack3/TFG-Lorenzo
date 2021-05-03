#' @title Computes a list of neighbouring antennas.
#'
#' @description Computes a list of pairs antennaID1-antennaID2 with neighbouring antennas. Two antennas are considered
#'   neighbours if their coverage areas overlap (i.e. their intersection is not void).
#'
#' @param coverarea a data.table object with two columns: 'antennaID' and 'cell'. The first column contains the ID of
#'   each antenna while the second one contains an \code{sp} object that represents the coverage area of the corresponding
#'   antenna. It is obtained by calling \code{readCells()} function.
#'
#' @return A data.table object with a single column called 'nei'. Each element of this column is a pair of the form
#'   antennaID1-antennaID2 where the two antennas are considered neighbours.
#'
#' @import data.table
#' @import rgeos
#' @export
antennaNeighbours <- function(coverarea) {
  
  antennas <- coverarea[['antennaID']]
  y <- data.table(expand.grid(antennas, antennas))
  y1 <-
    merge(
      y,
      coverarea,
      by.x = 'Var1',
      by.y = 'antennaID',
      all.x = TRUE,
      sort = TRUE
    )
  rm(y)
  y2 <-
    merge(
      y1,
      coverarea,
      by.x = 'Var2',
      by.y = 'antennaID',
      all.x = TRUE,
      sort = TRUE
    )
  rm(y1)
  colnames(y2) <- c('antennaID2', 'antennaID1', 'WKT1', 'WKT2')
  res <- vector(length = nrow(y2))
  for (i in 1:nrow(y2)) {
    p1 <- y2[[i, 'WKT1']]
    p2 <- y2[[i, 'WKT2']]
    res[i] <- gIntersects(p1, p2)
  }
  y2[, 'neighbour'] <- res
  y3 <- y2[neighbour == TRUE, .(antennaID1, antennaID2)]
  rm(y2)
  y3[, nei := do.call(paste, c(.SD, sep = "-"))]
  return(y3[, .(nei)])
}
