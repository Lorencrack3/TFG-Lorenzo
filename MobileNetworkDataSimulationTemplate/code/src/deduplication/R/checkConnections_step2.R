#' @title Make the second step of the checks between the connections and the emission model.
#'
#' @description This function obtains the results of the second step of checks
#' for the compatibility between the connections observed and the emission model. 
#' It checks that the transitions from one connections to the following is
#' coherent with the probabilities in the emission matrix. In those cases when 
#' the transition is not compatible, a solution to fit the model is given 
#' with the needed time padding coefficient.
#'
#' @param emissionProbs the matrix of emission probabilities.
#' 
#' @param connections a matrix with the connections, a row per device and
#' a column per time.
#' 
#' @param gridParams a list with the parameters of the grid: nrow, ncol, tileX, tileY.
#' 
#' @return A \code{list} with two elements, one is the information about the checks
#' with a vector containing all the time padding coefficients and the second 
#' is the matrix of connections with the time padding done by using the maximum
#' coefficient.
#' 
#' @examples
#' \dontrun{
#' }
#'
#' @import data.table 
#' 
#' @include distance_coverArea.R tileEquivalence.R
#' 
#' @export
checkConnections_step2 <- function(emissionProbs, connections, gridParams){
  
  # Check if emissionProbs meets the restrictions
  if (!all(round(apply(emissionProbs, 1, sum), 6)==1 || round(apply(emissionProbs, 1, sum), 6)==0)) {
    
    stop(paste0('[checkConnections_step2] There is some problem of consistency in emissionProbs\n'))
  }
  
  # Check if all observedValues are in emissionProbs
  columnsLack <- setdiff(connections, colnames(emissionProbs))
  columnsLack <- columnsLack[!is.na(columnsLack)]
  if (length(columnsLack) > 0) {
    
    stop(paste0('[checkConnections_step2] The following observedValues are missing in emissionProbs: ', 
                paste0(columnsLack, collapse = ', '),  '.\n'))
    
  }
  
  nrow_grid  <- gridParams$nrow
  ncol_grid  <- gridParams$ncol
  tile_sizeX <- gridParams$tileX
  tile_sizeY <- gridParams$tileY
  
  tileEquiv.dt <- data.table(tileEquivalence(ncol_grid, nrow_grid))
  
  centroidCoord_x  <- (0.5 + 0:(ncol_grid - 1)) * tile_sizeX
  centroidCoord_y  <- (0.5 + 0:(nrow_grid - 1)) * tile_sizeY
  centroidCoord.dt <- as.data.table(expand.grid(centroidCoord_y, centroidCoord_x))[
      , tile := 0:((ncol_grid * nrow_grid) - 1)][tileEquiv.dt, on = 'tile']
  setnames(centroidCoord.dt, c("Var1", "Var2"), c("centroidCoord_x", "centroidCoord_y"))
  
  coef_pad <- vector(length = nrow(connections), mode = "integer")
  
  for(i in 1:nrow(connections)){
    
    # cat(i)
    # cat(",\n")
    
    obs1 <- as.numeric(connections[i,])
    obs2 <- obs1[-1]
    aux.mt <- rbind(obs1[-length(obs1)], obs2)
    aux.mt <- rbind(aux.mt, aux.mt[1, ] - aux.mt[2, ])
    aux.mt <- aux.mt[, -which(aux.mt[3,] == 0), drop = FALSE]
    
    # For the direct transitions between antennas (without NA)
    pairs.dt <- unique(data.table(t(aux.mt[c(1:2), which(!is.na(aux.mt[3,]))])))
    setnames(pairs.dt, c("from", "to"))
    
    if(nrow(pairs.dt) > 0){
      
      # cat("check1:\n")
      # print(pairs.dt)
      pairs.dt[, dist_tiles := mapply(distance_coverArea,
                                      from, to, 
                                      MoreArgs = list(emissionProbs, centroidCoord.dt, tile_sizeX),
                                      SIMPLIFY = TRUE)]
      
      max_dist <- max(pairs.dt[, dist_tiles])
      if(max_dist > 0) coef_pad[i] <- max_dist + 1
      if(max_dist <= 0) coef_pad[i] <- 1
    }
    if(nrow(pairs.dt) == 0){
      coef_pad[i] <- 1
    }
    
    # For the transitions made with some NA 
    transNA.dt <- data.table(t(aux.mt[c(1:2), which(is.na(aux.mt[3,]))]))
    setnames(transNA.dt, c("from", "to"))
    if(nrow(transNA.dt) > 0){
      transNA.dt[!is.na(from), a:=1][!is.na(to), a:=2]
      # Drop the transitions from NA at the begining and to NA at the end
      aux_check <- transNA.dt[!is.na(a),a]
      if(aux_check[1] == 2){
        index_check1 <- transNA.dt[a == 1, which = TRUE]
        if(length(index_check1) > 0){
          keep <- index_check1[1]:nrow(transNA.dt)
        }
        if(length(index_check1) == 0) keep <- 0
        transNA.dt <- transNA.dt[keep]
      }
      if(aux_check[length(aux_check)] == 1 & nrow(transNA.dt) > 0){
        index_check2 <- transNA.dt[a == 2, which = TRUE]
        if(length(index_check2) > 0){
          keep <- 1:index_check2[length(index_check2)]
        }
        if(length(index_check2) == 0) keep <- 0
        transNA.dt <- transNA.dt[keep]
      }
      if(nrow(transNA.dt) > 0){
        pairsNA.dt <- data.table(from = transNA.dt[a==1, from], to = transNA.dt[a==2, to])
        NAtimes <- rle(is.na(transNA.dt[,from]))
        pairsNA.dt[, numNA := NAtimes$lengths[NAtimes$values]]
        # cat("check2:\n")
        # print(pairsNA.dt)
        pairsNA.dt[, dist_tiles := mapply(distance_coverArea,
                                          from, to, 
                                          MoreArgs = list(emissionProbs, centroidCoord.dt, tile_sizeX),
                                          SIMPLIFY = TRUE)]
        
        incre_coef <- max(pairsNA.dt[, dist_tiles] - pairsNA.dt[, numNA])
        if(incre_coef > 0){
          coef_pad[i] <- coef_pad[i] + (incre_coef + 1)
        }
      }
    }
    
  }
  
  pad_coef <- max(coef_pad)
  connections_pad <- matrix(NA, nrow = nrow(connections), 
                            ncol = pad_coef * ncol(connections))
  connections_pad[, seq(1, ncol(connections_pad), by = pad_coef)] <- connections
  
  output <- list(coef_pad = coef_pad, connections_pad = connections_pad)
  
  return(output)
  
}
  
  