#' @title Computes the probability distribution for Delta X and Delta Y.
#'
#' @description Computes the probability distribution for Delta X and Delta Y where Delta X and Delta Y are the
#'   differences between the centroids of the tiles on OX and OY. For a detailed explanation of what Delta X and Delta Y
#'   mean please reffer to
#'   \href{https://webgate.ec.europa.eu/fpfis/mwikis/essnetbigdata/images/f/fb/WPI_Deliverable_I3_A_proposed_production_framework_with_mobile_network_data_2020_05_31_draft.pdf}{WPI
#'    Deliverable 3} The posterior location probability for each tile and each device comes from \pkg{destim}.
#'
#' @param xp1 A data.table object with 3 columns representing the x and y coordinates of the centroids of tiles in the
#'   grid and the location probability for each corresponding tile for the first device.
#'
#' @param xp2 a data.table object with 3 columns representing the x and y coordinates of the centroids of tiles in the
#'   grid and the location probability for each corresponding tile for the second device.
#'
#' @return a list with two data.table objects, the first is the probability distribution for Delta X and the second for
#'   Delta Y. The elements are named DeltaX and DeltaY.
#'
#' @import data.table
#' @export
buildDeltaProb <-function(xp1, xp2) {
  
  colnames(xp1)<-c('x1', 'y1', 'p1')
  colnames(xp2)<-c('x2', 'y2', 'p2')
  
  deltaX <- CJ.table(xp1[,.(x1,p1)],xp2[,.(x2,p2)])
  deltaX <- deltaX[, delta:=x1-x2]
  deltaX <- deltaX[, p:=p1*p2]
  deltaX <- deltaX[, .(delta, p)]
  deltaX<-deltaX[, lapply(.SD,sum), by=.(delta)]
  
  deltaY <- CJ.table(xp1[,.(y1,p1)],xp2[,.(y2,p2)])
  deltaY <- deltaY[, delta:=y1-y2]
  deltaY <- deltaY[, p:=p1*p2]
  deltaY <- deltaY[, .(delta, p)]
  deltaY<-deltaY[, lapply(.SD,sum), by=.(delta)]
  
  return(list(deltaX=deltaX, deltaY=deltaY))
}


CJ.table <- function(X,Y) {
   setkey(X[,c(k=1,.SD)],k)[Y[,c(k=1,.SD)],allow.cartesian=TRUE][,k:=NULL]
}
