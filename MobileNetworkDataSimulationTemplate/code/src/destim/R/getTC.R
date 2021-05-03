# Gets the column of the transition matrix corresponding to
# a point from a rectangular grid according to the mask

# point : we are setting the transition probabilities
# starting from this point. It corresponds to the column
# (point[2] - 1 ) * size[1] + point[1] of the transition matrix.
#
# size : size of the grid. It is expected to be a 2 dimensional
# vector, containing the number of rows and columns repectively.
#
# mask : transition probabilities to the contiguous tiles. It
# is expected to be a 3x3 matrix where the (2,2) element
# represents the probability of staying in the current tile.

getTC <- function(point, size, mask) {
  output <- matrix(0, nrow = size[1], ncol = size[2])

  # Check if the point is in the vertices of the grid.
  if (all(point == 1))
    output[1:2,1:2] <- mask[2:3,2:3]
  else if (point[1] == 1 && point[2] == size[2])
    output[1:2, (size[2] - 1):size[2]] <- mask[2:3,1:2]
  else if (point[1] == size[1] && point[2] == 1)
    output[(size[1] - 1):size[1], 1:2] <- mask[1:2, 2:3]
  else if (point[1] == size[1] && point[2] == size[2])
    output[(size[1] - 1):size[1],
           (size[2] - 1):size[2]] <- mask[1:2, 1:2]

  # Check if the point is in the edges of the grid.
  else if (point[1] == 1)
    output[1:2, (point[2] - 1):(point[2] + 1)] <- mask[2:3,]
  else if (point[1] == size[1])
    output[(size[1] - 1):size[1],
           (point[2] - 1):(point[2] + 1)] <- mask[1:2,]
  else if (point[2] == 1)
    output[(point[1] - 1):(point[1] + 1), 1:2] <- mask[,2:3]
  else if (point[2] == size[2])
    output[(point[1] - 1):(point[1] + 1),
           (size[2] - 1):size[2]] <- mask[,1:2]
  # Otherwise is in the interior
  else output[(point[1] - 1):(point[1] + 1),
         (point[2] - 1):(point[2] + 1)] <- mask

  output <- as.vector(output)
  output <- output / sum(output)
  return(output)
}
