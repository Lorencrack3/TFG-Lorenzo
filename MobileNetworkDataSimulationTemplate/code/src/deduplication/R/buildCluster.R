#' @title Build a cluster of working nodes.
#'
#' @description Build a cluster of working nodes using the \code{parallel} package. On Unix-like operating systems it
#'   builds as "FORK" cluster while on Windows it builds a "SOCKET" cluster. Ther number of working nodes equals the
#'   number of available (logical) cores as it is detected by the call of \code{detectCores()} function. On Windows, it
#'   also exports a list of variables and an environment passed as arguments.
#'
#' @param varlist A list of variables to be exported to every working node.
#'
#' @param env An environment to be exported to every working node.
#' 
#' @return a cluster object.
#'
#' @import data.table
#' @import destim
#' @import parallel
#' @import doParallel
#'   
buildCluster <- function(varlist, env) {
  
  if (Sys.info()[['sysname']] == 'Linux' |
      Sys.info()[['sysname']] == 'Darwin') {
    cl <- makeCluster(detectCores(), type = "FORK")
  } else {
    cl <- makeCluster(detectCores())
    clusterEvalQ(cl, library("destim"))
    clusterEvalQ(cl, library("data.table"))
    clusterEvalQ(cl, library("Matrix"))
    clusterExport(cl, varlist, envir = env)
  }
  return (cl)
}
