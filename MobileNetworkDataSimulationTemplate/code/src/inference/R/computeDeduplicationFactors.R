#' @title  Computes the deduplication factors for each region
#'
#' @param dupFileName The name of the file with the duplicity probabilities. This file is the output of the
#'   \code{deduplication} package.
#'
#' @param regsFileName The name of the .csv file defining the regions. It has two columns: \code{ tile, region}. The
#'   first column contains the IDs of each tile in the grid while the second contains the number of a region. This file
#'   is defined by the user and it can be created with any text editor.
#'
#' @param postLocPath The path to the location where the posterior location probabilities are stored.
#'
#' @param postLocPrefix The file name prefix of the files with posterior location probabilities for each device. The
#'   whole file name is composed by a concatenation of \code{prefixName}, \code{_} and \code{deviceID}.The extension of
#'   these files is \code{.dt.csv}
#'
#'
#' @return A data.table object with the deduplication factors for each region.
#'
#'
#' @import data.table
#' @import parallel
#' @import deduplication
#' @export
computeDeduplicationFactors <- function(dupFileName, regsFileName, postLocPrefix, postLocPath) {
    
    if (!file.exists(dupFileName))
        stop(paste0(dupFileName, " does not exists!"))
    
    duplicity <- fread(
        dupFileName,
        sep = ',',
        header = TRUE,
        stringsAsFactors = FALSE
    )
    setnames(duplicity, c('device', 'dupP'))
    
    if (!file.exists(regsFileName))
        stop(paste0(regsFileName, " does not exists!"))
    
    regions <- fread(
        regsFileName,
        sep = ',',
        header = TRUE,
        stringsAsFactors = FALSE
    )
    
    devices<-sort(duplicity$device)
    

    ndevices<-length(devices)
    
    postLocPath <- postLocPath
    postLocPrefix <- postLocPrefix
    
    cl <- buildCluster( c('postLocPath', 'postLocPrefix', 'devices') , env=environment() )
    ich<-clusterSplit(cl, 1:ndevices)
    res<-clusterApplyLB(cl, ich, doRead, postLocPath, postLocPrefix, devices)
    postLocProb0 <- rbindlist((res))
    stopCluster(cl)
    setnames(postLocProb0,'V2', 'device')

    postLocProb0 <- merge(postLocProb0, regions, by = 'tile', all.x = TRUE)[
            , list(prob = sum(probL)), by = c('device', 'region')]

    
    omega_region <- merge(postLocProb0, duplicity, by='device')[
        , list(omega1 = sum(prob * (1 - dupP)) / sum(prob), 
               omega2 = sum(prob * dupP) / sum(prob)), by = 'region'][
        order(region)]
    
    return (omega_region)
    
}

doRead <- function(ichunks, postLocPath, postLocPrefix, devices) {
    local_postLoc <- NULL
    for( i in ichunks) {
        tmp <- readPostLocProb(postLocPath, postLocPrefix, devices[i])
        tmp <- tmp[time==0][,time:=NULL][,device := NULL]
        nr <- nrow(tmp)
        tmp <- cbind(tmp, rep(devices[i], times = nr))
        local_postLoc <-rbind(local_postLoc, tmp)
    }
    return (local_postLoc)
}
