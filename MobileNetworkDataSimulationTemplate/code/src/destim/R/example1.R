#' This is just an example
#'
#' @examples
#' # Local probabilities of transition
#' mask <- matrix(c(0.00001, 0.00320, 0.00001,
#'                  0.00320, 0.98716, 0.00320,
#'                  0.00001, 0.00320, 0.00001), ncol = 3)
#' #Complete transition matrix
#' P <- createTM(c(20,20), mask)
#' #Towers position
#' data(towers)
#' # Function S for Tennekes's model
#' S <- function(x) if (x > 5) return(0) else return(20*log(5/x))
#' # Complete events matrix
#' E <- createEM(c(20,20), towers, S)
#' # Load detection events
#' data(events)
#' # Estimate filtered states
#' FS <- fstates(P, E, events)
#' # Estimate smooth states
#' SS <- sstates(P, E, events, FS)
#' # Load some required packages
#' library(raster)
#' library(RColorBrewer)
#' library(animation)
#' library(ggplot2)
#' # Create the animations
#' # Filtered states
#' saveGIF({
#' pal <- colorRampPalette(c("white","red"))
#' ani.options(interval = 0.02)
#' for (i in 1:198) {
#'  plot(raster(matrix(FS[,i], ncol = 20)), zlim = c(0,1), col = pal(100))
#'  ani.pause()
#' }},  movie.name = 'filteredest.gif')
#' # Smooth states
#' saveGIF({
#' pal <- colorRampPalette(c("white","red"))
#' ani.options(interval = 0.02)
#' for (i in 1:198) {
#'  plot(raster(matrix(SS[,i], ncol = 20)), zlim = c(0,1), col = pal(100))
#'  ani.pause()
#' }},  movie.name = 'smoothest.gif')
#' # The matrix GRID relates the states with coordinates
#' GRID <- matrix(1, nrow = 2, ncol = 400)
#' GRID[1,] <- rep(1:20,20)
#' GRID[2,] <- rep(1:20, each = 20)
#' # Calculate square distance mean
#' fdist <- sapply (1:198, function (x)
#'           sum(apply(GRID - matrix(c(x %/% 10 + 1,14), nrow = 2, ncol = 400), 2,norm,type = "2") * FS[,x]))
#' sdist <- sapply (1:198, function (x)
#'           sum(apply(GRID - matrix(c(x %/% 10 + 1,14), nrow = 2, ncol = 400),2,norm,type = "2") * SS[,x]))
#' dists <- data.frame(T = 1:198, fdist = fdist, sdist = sdist)
#' ggplot() +
#' geom_line(data = dists, aes(x = T, y = sdist), colour = "green") +
#' geom_line(data = dists, aes(x = T, y = fdist), colour = "red")
#'
example1 <- function() {}
