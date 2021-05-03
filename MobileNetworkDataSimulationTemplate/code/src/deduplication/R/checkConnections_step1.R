#' @title Make the first step of the checks between the connections and the emission model.
#'
#' @description This function obtains the results of the first step of checks
#' for the compatibility between the connections observed and the emission model. 
#' It checks if each individual connection is coherent with the probabilities 
#' in the emission matrix. In those cases when that connection is not compatible, 
#' the connection is imputed with a missing value.
#'
#' @param connections a matrix with the connections, a row per device and
#' a column per time.
#' 
#' @param emissionProbs a matrix of emission probabilities resulting from 
#' the function getEmissionProbs in deduplication package.
#' 
#' @return A \code{list} with two elements, one is the information about the checks
#' and the second is the matrix of connections with imputation in the observations
#' not compatible with the model in emissionProbs.
#'
#' @examples
#' \dontrun{
#' }
#'
#' @import data.table
#'
#' @export
checkConnections_step1 <- function(connections, emissionProbs){
  
  
  checkE <- apply(emissionProbs, 2, sum)
  
  infoCheck1.dt <- data.table()
  connections2 <- connections
  
  for(i in 1:nrow(connections)){
    
    events_i <- unlist(connections[i,])
    condition_i <- checkE[events_i]==0
    condition_i[is.na(condition_i)] <- TRUE
    nocon_i <- events_i[condition_i]
    time_i <- c(1:ncol(connections))[condition_i]
    
    if(length(nocon_i) > 0){
      aux.dt <- data.table(deviceID = devices[i], 
                           antennaID = nocon_i, 
                           time = time_i)
      infoCheck1.dt <- rbindlist(list(infoCheck1.dt, aux.dt))[!is.na(antennaID)]
      connections2[i, condition_i] <- NA
    }
    
  }
  
  output <- list(infoCheck_step1 = infoCheck1.dt, connectionsImp = connections2)
  return(output)
}


