#' @title  Builds a concatenation between centroids and location probabilities.
#'
#' @description
#'
#' @param centroids The centroids of the tiles in the grid. The centroids are computed by calling
#'   \code{buildCentroids()} function.
#'
#' @param postLoc The posterior location probabilities for each tile. 
#'
#' @param t The time instant for which the function uses the posterior location probabilities.
#'
#' @return A data.table object with 3 columns: x and y coordinates of the centroid of a tile and the location
#'   probability for that tile at a specific time instant.
#'
#' @export
buildCentroidProbs <-function(centroids, postLoc, t) {
  
  xp<-list(length = length(postLoc))
  for(i in 1:length(postLoc)) {
    tmp<-cbind(centroids[,.(x,y)], postLoc[[i]][,t])
    #coln<-colnames(tmp)
    tmp<-tmp[V2!=0][,]
    xp[[i]]<-tmp
  }

  return (xp)
  
}