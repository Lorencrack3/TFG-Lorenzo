# Creates the events matrix for a rectangular grid
# according to the location of the towers and the S function

createEM <- function(size, towers, S) {
  output <- matrix(0, nrow = size[1] * size[2],
                   ncol = ncol(towers))

  for (i in 1:size[1])
    for (j in 1:size[2])
      output[(j - 1) * size[2] + i, ] <-
        getER(c(i,j), towers, S)

  return(output)
}
