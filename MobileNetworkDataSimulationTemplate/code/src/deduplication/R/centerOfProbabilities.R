#' @title Computes the center of probabilities (as center of mass) for each posterior probability cloud, at each time
#'   instant, for a device.
#'
#' @description Computes the center of probabilities (as center of mass) for each posterior probability cloud, at each
#'   time instant, for a device. It uses the results obtained from \pkg{destim} package.
#'
#' @param centroid A data.table object with the centroids for each tile in the grid.
#'
#' @param postLocProb a Matrix object with the posterior location probabilities computed by \code{destim} for a device.
#'   The number of rows equals the number of tiles and the number of columns equals the number of successive time
#'   instants. An element \code{postLocProb[i,j]} represents the posterior location probability corresponding to tile
#'   \code{i} at the time instant \code{j}.
#'
#' @return A table with two columns: centerProb_x and centerProb_y representing the x and y coordinates of the centers
#'   of probabilities.
#'
#' @import data.table
#' @export
centerOfProbabilities <- function(centroid, postLocProb) {
  
  nc <- nrow(centroid)
  if(length(postLocProb) == 1) {
    if(postLocProb * nc != 1)
      stop('in case all tiles have the same posterior location probability, the sum of all these values must be 1')

  }
  else { 
    nw <- nrow(postLocProb)
    if (nw != 1 & nw != nc)
      stop('postLocProb must have values for all tiles.')
  }

  cp_x <- vector(mode ='numeric', length = ncol(postLocProb))
  cp_y <- vector(mode = 'numeric', length = ncol(postLocProb))
  for (i in 1:ncol(postLocProb)) {
    cp_x[i] <- sum(postLocProb[, i] * centroid[, 1])# / sum(postLocProb[, i]) #inutila impartirea (numitorul e 1)
    cp_y[i] <- sum(postLocProb[, i] * centroid[, 2])# / sum(postLocProb[, i]) #inutila impartirea (numitorul e 1)
  }
  cp <- cbind(cp_x, cp_y)
  colnames(cp) <- c('centerProb_x', 'centerProb_y')
  return(cp)
  
}