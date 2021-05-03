#' Creates a basic rectangle HMM model
#'
HMMrectangle <- function(x,y) {
  x <- as.integer(x)
  y <- as.integer(y)
  if ((x <= 0) || (y <= 0))
    stop("Both dimensions have to be integer positive numbers.")
  TL <- matrix(0, nrow = 2, ncol = 9*x*y - 6*x - 6*y + 4)
  k <- 1
  for (i in 1:x)
    for (j in 1:y)
      for (newi in (max(i - 1, 1)):(min(i + 1, x)))
        for (newj in (max(j - 1, 1)):(min(j + 1, y))) {
            TL[, k] <- c((j - 1) * x + i, (newj - 1) * x + newi)
            k <- k + 1
        }
  newmodel <- HMM(x*y, TL)
  TL <- TL[2,] - TL[1,]
  # Dirty code to deal with x == 2
  if (x == 2) {
    if (y == 2)
      TL[which(apply(transitions(newmodel), 2, function(x)
         all (x == c(2L,3L)) | all (x == c(3L,2L))))] <- c(5L, -5L)
    if (y > 2) {
      TL[which(apply(transitions(newmodel), 2, function(x)
        (x[2] - x[1] == 1) & (x[1] %% 2 == 0)))] <- 5L
      TL[which(apply(transitions(newmodel), 2, function(x)
        (x[2] - x[1] == -1) & (x[1] %% 2 == 1)))] <- -5L
    }
  }
  if (min(x,y) == 1)
    CT <- matrix(0, nrow = length(TL) - length(unique(TL)) -
                   (max(x,y) > 2),
                 ncol = ntransitions(newmodel) + 1)
  else if (max(x,y) == 2)
    CT <- matrix(0, nrow = 7,
               ncol = ntransitions(newmodel) + 1)
  else
    CT <- matrix(0, nrow = length(TL) - length(unique(TL)) - 1 -
                             (min(x, y) > 2),
                 ncol = ntransitions(newmodel) + 1)
  k <- 1
  for (i in unique(TL))
    if (i != 0) {
      j <- which(TL == i)
      if (length(j) > 1)
        for (l in 1:(length(j) - 1)) {
          CT[k, j[l]] <- 1
          CT[k, j[l + 1]] <- -1
          k <- k + 1
        }
    }
  stateclasses <- sapply(1:nstates(newmodel), function(n)
                                sum(transitions(newmodel)[1,] == n))
  for (j in unique(stateclasses)) {
  previous <- 0
  for (i in which(stateclasses == j))
    if (previous == 0)
      previous <- i
    else {
      CT[k, c((transitions(newmodel)[1, ] ==
                 previous) & (TL == 0), FALSE)] <- 1
      CT[k, c((transitions(newmodel)[1, ] ==
                 i) & (TL == 0), FALSE)] <- -1
      k <- k + 1
    }
  }
  newmodel <- addconstraint(newmodel, CT)
  hvmov <- which(TL %in% c(1, -1, x, -x))[1]
  for (i in setdiff(c(1, -1, x, -x), TL[hvmov]))
    if (any(TL == i))
      newmodel <- addconstraint(newmodel, c(hvmov, which(TL == i)[1]))
  diagmov <- which(!(TL %in% c(1,-1, 0, x, -x)))[1]
  if (!is.na(diagmov))
    for (i in setdiff(unique(TL), c(1,-1, 0, x, -x, TL[diagmov])))
      newmodel <- addconstraint(newmodel,
                                c(diagmov, which(TL == i)[1]))
  newmodel$constraints <-
    Reduce(rbind, simplifyconstraints (newmodel$constraints))
  return(newmodel)
}
