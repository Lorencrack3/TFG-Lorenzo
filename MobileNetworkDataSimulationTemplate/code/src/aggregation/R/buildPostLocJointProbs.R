#'@title Builds the table with joint posterior location probabilities from
#'  simulated data.
#'
#'@description This is a utility function with the purpose to build the files
#'  with joint posterior location probabilities for each device using the
#'  \code{destim} package an a series of additional data. This function was used
#'  to build the example data set included in this package The raw data were
#'  produced using the simulation software. The joint posterior location
#'  probabilities file names starts with \code{postLocJointProbDevice_} followed for
#'  the \code{deviceID} and the extension \code{.dt.csv}.
#'
#'@param path The path of the directory where the outfiles are saved.
#'
#'@param simFileName The name of the file with the general simulation parameters.
#'  parameters.
#'
#'@param gridFileName The name of the file with the grid parameters.
#'
#'@param eventsFileName The name of the file containing the networks events.
#'
#'@param signalFileName The name of the file with the signal strength/quality.
#'
#'
#'@import data.table
#'@import destim
#'@import deduplication
#'@export
buildPostLocJointProbs <-function(path, simFileName, gridFileName, eventsFileName, signalFileName) {
    if(!file.exists(simFileName))
        stop(paste0(simFileName, " doesn't exists"))
    
    if(!file.exists(gridFileName))
        stop(paste0(gridFileName, " doesn't exists"))
    
    if(!file.exists(eventsFileName))
        stop(paste0(eventsFileName, " doesn't exists"))
    
    if(!file.exists(signalFileName))
        stop(paste0(signalFileName, " doesn't exists"))
    
    simParams <- readSimulationParams(simFileName)
    gridParams <- readGridParams(gridFileName)
    events <- readEvents(eventsFileName)
    devices <- getDeviceIDs(events)
    connections <- getConnections(events)
    emissionProbs <- getEmissionProbs(gridParams$nrow, gridParams$ncol, signalFileName, simParams$conn_threshold)
    model <- getGenericModel(gridParams$nrow, gridParams$ncol, emissionProbs)
    eq<-as.data.table(tileEquivalence(gridParams$nrow, gridParams$ncol))
    for (j in 1:length(devices) ){
        modeli <- fit(model, connections[j,], init = TRUE, method = "solnp")
        B <- scpstates(modeli, connections[j,])
        out<-transform_postLocJointProb(postJointLocP = B, observedValues=connections[j,], 
                                        times=seq(from=simParams$start_time, to=simParams$end_time, by=simParams$time_increment), 
                                        t_increment = simParams$time_increment, 
                                        ntiles=gridParams$nrow * gridParams$ncol, 
                                        tileEquiv.dt = eq, devID=devices[j], 
                                        sparse_postJointLocP.dt = TRUE)
        out[,rasterCell_from:=NULL]
        out[,rasterCell_to:=NULL]
        out[,device:=NULL]
        write.table(out, file=paste0(path, "/postLocJointProbDevice_", devices[j],".dt.csv"), sep=",", row.names = FALSE)
    }
}


transform_postLocJointProb <- function(postJointLocP, 
                                       observedValues,
                                       times, 
                                       t_increment, 
                                       ntiles, 
                                       pad_coef = 1, 
                                       tileEquiv.dt, 
                                       devID,
                                       sparse_postJointLocP.dt = FALSE){
    
    
    if(!is.null(postJointLocP)){
        
        postJointLocP0 <- copy(postJointLocP)
        postJointLocP0 <- as(postJointLocP0, "dgTMatrix")
        postJointLocP.dt <- data.table(rasterCell_from = postJointLocP0@i+1, 
                                       j = postJointLocP0@j+1, 
                                       probL = postJointLocP0@x)
        postJointLocP.dt <- postJointLocP.dt[, time_from := floor((j - 1) / ntiles) + 1][, rasterCell_to := j - ((time_from - 1) * ntiles)]
        postJointLocP.dt <- postJointLocP.dt[, j := NULL]
        
        postJointLocP.dt <- postJointLocP.dt[time_from %in% seq(1, length(observedValues), by = pad_coef)]
        names(times) <- sort(unique(postJointLocP.dt$time_from))
        postJointLocP.dt[, time_from := times[as.character(time_from)]][
            , device := devID]
        postJointLocP.dt <- merge(postJointLocP.dt, tileEquiv.dt, by.x = 'rasterCell_from', by.y = 'rasterCell')
        setnames(postJointLocP.dt, 'tile', 'tile_from')
        postJointLocP.dt <- merge(postJointLocP.dt, tileEquiv.dt, by.x = 'rasterCell_to', by.y = 'rasterCell')
        setnames(postJointLocP.dt, 'tile', 'tile_to')
        
        if(!sparse_postJointLocP.dt){
            allTiles.dt <- data.table(expand.grid(c(1:ntiles), c(1:ntiles)))
            setnames(allTiles.dt, c("Var1", "Var2"), c("rasterCell_from", "rasterCell_to"))
            allTiles.dt[, aux := 1]
            times.dt <- data.table(time_from = sort(unique(postJointLocP.dt$time_from)), aux = 1)
            allTilesTimes <- merge(allTiles.dt, times.dt, by = "aux", allow.cartesian = TRUE)
            allTilesTimes[, aux := NULL][, device := devID]
            allTilesTimes <- merge(allTilesTimes, tileEquiv.dt, by.x = 'rasterCell_from', by.y = 'rasterCell')
            setnames(allTilesTimes, 'tile', 'tile_from')
            allTilesTimes <- merge(allTilesTimes, tileEquiv.dt, by.x = 'rasterCell_to', by.y = 'rasterCell')
            setnames(allTilesTimes, 'tile', 'tile_to')
            
            postJointLocP.dt <- merge(allTilesTimes, postJointLocP.dt, 
                                      all.x = TRUE, by = intersect(names(allTilesTimes), names(postJointLocP.dt)))
            rm(allTilesTimes)
            postJointLocP.dt[is.na(probL), probL := 0]
            
        }
        
        postJointLocP.dt <- postJointLocP.dt[, time_to := time_from + t_increment]
        postJointLocP.dt[probL<0, probL := 0]
        setcolorder(postJointLocP.dt, c('time_from', 'time_to', 'tile_from', 'tile_to', 
                                        'probL', 'device', 'rasterCell_from', 'rasterCell_to'))
        
    }
    
    return(postJointLocP.dt)
    
}
