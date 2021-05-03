#' @title Reads the coverage areas of antennas.
#'
#' @description Reads the coverage areas of antennas from a .csv file.
#'
#' @param cellsFileName It is the name of the file where the coverage areas of antennas are to be found. The data have two
#'   columns, the first one is the antenna ID and the second one is a WKT string representing a polygon (i.e. it should
#'   start with the word POLYGON) which is the coverage area of the corresponding antenna. This area is also called the
#'   antenna cell.
#'
#' @param simulatedData If TRUE it means that the file with the coverage areas is produced by the data simulator
#'
#' @return A data.table object with 2 columns: the antenna ID and an \code{sp} geometry object which is the coverage
#'   area of the corresponding antenna.
#'
#' @import data.table
#' @import rgeos
#' @export
readCells <- function (cellsFileName, simulatedData = TRUE) {
  if (!file.exists(cellsFileName))
    stop(paste0(cellsFileName, " does not exist!"))
  
  if (simulatedData) {
    coverArea <-
      fread(
        cellsFileName,
        sep = '\n',
        header = TRUE,
        stringsAsFactors = FALSE
      )
    setnames(coverArea, 'lines')
    coverArea[, antennaID := tstrsplit(lines, split = ',POLYGON')[[1]]]
    coverArea[, wkt := substring(lines, regexpr("POLYGON", lines))]
    coverArea <- coverArea[, c('antennaID', 'wkt'), with = FALSE]
    coverArea[, 'cell'] <- sapply(coverArea[['wkt']], function(wkt) {
      polygon <- readWKT(wkt)
      return(polygon)
    })
    return (coverArea[, -2])
  }
  else {
    print("Read real mobile network cell file not implemented yet!")
    return (NULL)
  }
}
