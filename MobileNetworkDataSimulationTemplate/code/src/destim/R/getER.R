# Gets the row of the events matrix corresponding to
# a point from a rectangular grid according to the location
# of the towers and the S function

getER <- function(point, towers, S) {
  towers <- towers - matrix(point, nrow = 2, ncol = ncol(towers))
  distances <- apply(towers, 2, norm, type = "2")
  output <- sapply(distances, S)
  output <- output / sum(output)
  return(output)
}
