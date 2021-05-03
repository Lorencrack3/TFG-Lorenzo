# Creates the transition matrix of a rectangular grid
# according to the given mask

createTM <- function(size, mask) {
  output <- matrix(0, nrow = size[1] * size[2],
                   ncol = size[1] * size[2])
  if (length(dim(mask)) == 2)
    for (i in 1:size[1])
      for (j in 1:size[2]) {
        output[, (j - 1) * size[2] + i] <-
          getTC(c(i,j), size, mask)
      } else for (i in 1:size[1])
      for (j in 1:size[2]) {
        output[, (j - 1) * size[2] + i] <-
          getTC(c(i,j), size, mask[i,j,,])
      }
  return(output)
}
