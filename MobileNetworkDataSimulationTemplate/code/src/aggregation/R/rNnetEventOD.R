#'@title Generates random value according to a Poisson multinomial distribution
#'  needed to estimate the origin destination matrices.
#'
#'@description Generates random value according to a Poisson multinomial
#'  distribution needed to estimate the origin destination matrices. This is a
#'  high level function, the only one to be called by users to estimate the
#'  number of individuals going from one region to another. The actual
#'  computations are performed using a parallelization (transparent to the
#'  users) which uses the whole number of (logical) cores.
#'
#'@param n The number of random values to be generated.
#'
#'@param dupFileName The name of the .csv file with the duplicity probability
#'  for each device. This is an output of the \code{deduplication} package.
#'
#'@param regsFileName The name of the .csv file defining the regions. It has two
#'  columns: \code{ tile, region}. The first column contains the IDs of each
#'  tile in the grid while the second contains the number of a region. This file
#'  is defined by the user and it can be created with any text editor.
#'
#'@param postLocJointPath The path where the files with the posterior location
#'  probabilities for each device can be found. A file with the location
#'  probabilities should have the name \code{prefix_ID.csv} where \code{ID} is
#'  replaced with the device ID and \code{prefix} is given as a parameter to
#'  this function.
#'
#'@param prefix A prefix that is used to compose the file name with posterior
#'  location probabilities.
#'
#'@return  A data table object with the following columns: \code{time_from,
#'  time_to, region_from, region_to, Nnet, iter}. The number of detected
#'  individuals moving from a region to another between two succesive time
#'  instants is given in column \code{Nnet} while the last column gives the
#'  index of the random value generated for this number.
#'
#'
#'
#'@import extraDistr
#'@import data.table
#'@import deduplication
#'@import parallel
#'
#' @include rNnetJoint_Event.R
#'
#'@export
rNnetEventOD <- function(n, dupFileName, regsFileName, postLocJointPath, prefix, seed = 123) {
    
    if (!file.exists(dupFileName))
        stop(paste0(dupFileName, " does not exists!"))
    
    dupProbs <- fread(
        dupFileName,
        sep = ',',
        header = FALSE,
        stringsAsFactors = FALSE
    )
    setnames(dupProbs, c('device', 'dupP'))
    devices<-as.numeric(dupProbs[,1][[1]])

        
    if (!file.exists(regsFileName))
        stop(paste0(regsFileName, " does not exists!"))
    
    regions <- fread(
        regsFileName,
        sep = ',',
        header = TRUE,
        stringsAsFactors = FALSE
    )
    
    ndevices <- nrow(dupProbs)
    postLocJoint<-NULL
    
    cl <- buildCluster(c('postLocJointPath', 'prefix', 'devices', 'readPostLocProb') , env=environment())
    clusterSetRNGStream(cl, iseed=seed)
    ichunks <- clusterSplit(cl, 1:ndevices)
    res <-
        clusterApplyLB(
            cl,
            ichunks,
            doRead,
            postLocJointPath,
            prefix,
            devices
        )
    for(i in 1:length(res))
        res[[i]]<-as.data.table(as.matrix(res[[i]]))
    result <- rbindlist(res)
    
    postLocJoint <- as.data.table(as.matrix(result))[
        , c('device', 'time_from', "time_to", "tile_from", "tile_to", "probL")]
#    setnames(postLocJoint, c('device', 'time_from', "time_to", "tile_from", "tile_to", "eventLoc"))

    times<-sort(unlist(unique(postLocJoint[,2])))
    time_increment <- times[2]-times[1]
    #rm(times)
    
    postLocJointReg <- merge(
        postLocJoint,
        regions, 
        by.x = 'tile_from', by.y = 'tile')
    setnames(postLocJointReg, 'region', 'region_from')

    postLocJointReg <- merge(
        postLocJointReg,
        regions,
        by.x = 'tile_to', by.y = 'tile')
    setnames(postLocJointReg, 'region', 'region_to')
    
    postLocJointReg <- postLocJointReg[, list(prob = sum(probL)), by = .(device, time_from, time_to, region_from, region_to)]
    
    dedupProbs2_1_Reg <- merge(
        postLocJointReg, dupProbs[, device := as.numeric(device)], 
        by = 'device', all.x = TRUE)
    
    dedupProbs1_1_Reg <- copy(dedupProbs2_1_Reg)[
        , singleP := 1 - dupP][
            , dupP := NULL]
    
    dedupProbs2_1_Reg[
        , prob := prob * dupP][
            , devCount := 0.5][
                , dupP := NULL]
    
    dedupProbs1_1_Reg[
        , prob := prob * singleP][
            , devCount := 1][
                , singleP := NULL]
    
    dedupProbsReg <- rbindlist(
        list(dedupProbs1_1_Reg, dedupProbs2_1_Reg))
    rm(dedupProbs1_1_Reg, dedupProbs2_1_Reg)

    time_from <- unique(unlist(dedupProbsReg[,time_from]))
    ichunks <- clusterSplit(cl, time_from)

    cellNames<-unname(sort(unlist(unique(regions[,2]))))

    n<-n
    #setnames(dedupProbsReg, c('region_from', 'region_to'), c('cell_from', 'cell_to'))

    clusterExport(cl, c('cellNames', 'n', 'time_increment', 'dedupProbsReg'), envir = environment())
    res <-
        clusterApplyLB(
            cl,
            ichunks,
            doOD2,
            n,
            cellNames,
            time_increment,
            dedupProbsReg
        )
    stopCluster(cl)
    result <- rbindlist(res)
    return(result)
}

doOD2 <- function(ichunks, n, cellNames, time_increment, dedupProbs){
    
    NnetReg <- dedupProbs[time_from %in% ichunks][
        , aggregation::rNnetJoint_Event(.SD, cellNames = cellNames, n=n ), by = 'time_from', .SDcols = names(dedupProbs)][
            , time_to := time_from + time_increment][
                , region_to := as.integer(region_to)][
                    , region_from := as.integer(region_from)]
    setcolorder(NnetReg, c('time_from', 'time_to', 'region_from', 'region_to', 'Nnet', 'iter'))
    return(NnetReg)
}

doRead <- function(ichunks, path, prefix, devices) {
    postLoc <- NULL
    for(i in ichunks) {
        l <- readPostLocJointProb(path, prefix, devices[i])
        #l <- cbind(rep(devices[i], times= nrow(l)), l )
        postLoc<-rbind(postLoc,l)
    }
    return(postLoc)    
}