#' @title Utility function used throughout the package.
#'
#' @description Contains a function to compute the Mode of a data set 
#'
#' @export
Mode = function(v){
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

