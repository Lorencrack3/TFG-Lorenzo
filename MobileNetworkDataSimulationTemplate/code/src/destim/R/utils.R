# These functions are for internal use and won't be exported in the future

funique <- function(x) {
  return(cppfunique(x, sqrt(.Machine$double.eps)))
}

funiqueind <- function(x) {
  return(cppfuniqueind(x, sqrt(.Machine$double.eps)))
}

simplifyconstraints <- function(CT) {
  trmatrix <- matrix(as.logical(diag(ncol(CT) - 1)),
                     ncol = ncol(CT) -1)
  eqcon <- apply(CT, 1, function(irow)
    (irow[ncol(CT)] == 0) &&
      (sum(irow != 0) == 2) &&
      (sum(irow) == 0))
  EQCT <- CT[eqcon, 1:(ncol(CT) - 1), drop = FALSE]
  trmatrixcols <- ncol(trmatrix)
  for (i in seq_len(nrow(EQCT))) {
    variables <- match(c(1, -1), EQCT[i, ], incomparables = 0)
    columns <- c(match(TRUE, trmatrix[variables[1], ]),
                 match(TRUE, trmatrix[variables[2], ]))
    columns <- sort(columns)
    trmatrix[,columns[1]] <- trmatrix[,columns[1]] |
      trmatrix[,columns[2]]
    trmatrix[, columns[2]] <-
      trmatrix[, trmatrixcols, drop = FALSE]
    trmatrixcols <- trmatrixcols - 1
  }
  trmatrix <- trmatrix[,1:trmatrixcols, drop = FALSE]
  trmatrix <- matrix(as.numeric(trmatrix), ncol = ncol(trmatrix))
  MCT <- CT[!eqcon, , drop = FALSE]
  if (nrow(MCT) > 0) {
    MCT <- cbind(MCT[, -ncol(MCT), drop = FALSE] %*% trmatrix,
               MCT[,ncol(MCT), drop = FALSE])
    MCT <- CT[!eqcon, , drop = FALSE][funiqueind(MCT), ]
  }
  EQCT <- cbind(EQCT, rep(0, nrow(EQCT)))
  return(list(EQCT,MCT))
}
