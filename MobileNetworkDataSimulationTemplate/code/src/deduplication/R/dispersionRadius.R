#' @title Computes the radius of dispersion for the probability cloud of a device at each time instant.
#'
#' @description  Computes the radius of dispersion for the probability cloud of a device at each time instant. The
#'   posterior location probabilities for a specific device, for each tile and each time instant, is provided as a
#'   parameter. The definition of the dispersion radius is given in
#'   \href{https://webgate.ec.europa.eu/fpfis/mwikis/essnetbigdata/images/f/fb/WPI_Deliverable_I3_A_proposed_production_framework_with_mobile_network_data_2020_05_31_draft.pdf}{WPI
#'    Deliverable 3}.
#'
#' @param centroid The centroids of all tiles in the grid.
#'
#' @param postLocProb The posterior location probabilities for a device, for each tile and each time instant.
#'
#' @param center The (mass) center of posterior location probabilities. If it is not provided, it is computed here.
#'
#' @param method The distance function to be used. The euclidean distance is the default method.
#'
#' @param p The power of the Minkowski distance. The default value is 2 (euclidean distance).
#'
#' @return A vector with dispersion radius for each time instant, for the device under consideration.
#'
#' @import data.table
#' @include centerOfProbabilities.R
#' @export
dispersionRadius <- function(centroid, postLocProb, center = NULL, method = 'euclidean', p = 2) {

    nc <- nrow(centroid)
    nl <- length(postLocProb)

    if(nl == 1) {
      if(postLocProb * nc != 1)
        stop('in case all tiles have the same posterior location probability, the sum of all these values must be 1')
    } else { 
      nw <- nrow(postLocProb)
      if (nw != 1 & nw != nc)
        stop('postLocProb must have values for all tiles.')
    }
    
    if (is.null(center)) {
      center <- centerOfProbabilities(centroid, postLocProb)
    }

    ntimes<-nrow(center)
    rd<-vector(length = ntimes)
    tmp <- as.matrix(centroid[,1:2])
    for(i in 1:ntimes) {
      coordMatrix <- rbind(center[i,], tmp)
      distance <- dist(coordMatrix, method = method, p = p)[1:nc]
      rd[i] <- sqrt(sum(postLocProb[,i] * distance**2))
    }
    return(rd)
}
