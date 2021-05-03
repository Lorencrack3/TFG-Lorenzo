####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
path_root         <- 'G:/GRUPO_BIGDATA/Proyecto_ESSNet Big Data II/Simulations/template'
path_source       <- file.path(path_root, 'code/src')
path_simConfig    <- file.path(path_root, 'data/simulatorConfig')
path_events       <- file.path(path_root, 'data/networkEvents')
path_eventLoc     <- file.path(path_root, 'data/eventLocProb')
path_resources    <- file.path(path_root, 'param/resources')
path_processParam <- file.path(path_root, 'param/process')
path_postLoc      <- file.path(path_root, 'data/postLocProb')
path_dupProb      <- file.path(path_root, 'data/dupProb')


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)
library(xml2)                 # to read simulatorConfig and parameters
library(tidyr)                # transform xml documents to tables (tibbles)
library(rgeos)                # spatial data
library(destim)
library(stringr)
library(ggplot2)              # graphics
library(pROC)                 # for ROC curve for duplicity test
library(deduplication)

# Function get_simConfig to read the input files of the simulator
source(file.path(path_source, 'get_simConfig.R'))
# Function get_simScenario to read the output files of the simulator
source(file.path(path_source, 'get_simScenario.R'))
# Function tileEquivalence to compute the equivalence between rastercell (R) and tiles (simulator)
source(file.path(path_source, 'tileEquivalence.R'))


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           LOAD PARAMETERS                              ####
parameters.xml <- as_list(read_xml(file.path(path_processParam, "parameters_process.xml")))

emission_model    <- parameters.xml$process_parameters$geolocation$emission_model[[1]]
geolocation_prior <- parameters.xml$process_parameters$geolocation$prior[[1]]
method            <- parameters.xml$process_parameters$deduplication$method[[1]]
lambda            <- eval(parse(
  text = parameters.xml$process_parameters$deduplication$lambda[[1]]))

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                              COMPUTATION                               ####

# 0. Read simulation params                                                 ####
simulationFileName <- file.path(path_simConfig, 'simulation.xml')
simParams          <- readSimulationParams(simulationFileName)

# 1. Read grid parameters                                                   ####
gridFileName <- file.path(path_resources, 'grid.csv')
gridParams   <- readGridParams(gridFileName)

# 2.Read network events                                                     ####
eventsFileName <- file.path(path_events, 'AntennaInfo_MNO_MNO1.csv')
events         <- readEvents(eventsFileName)

# 3. Set signal file name                                                   ####
### i.e. the file where the signal strength/quality for each tile in the grid is stored
signalFileName <- file.path(path_resources, 'SignalMeasure_MNO1.csv')

# 4. Set antenna file name                                                  ####
### , i.e. the file where the signal strength/quality for each tile in the grid is stored
antennasFileName <- file.path(path_simConfig, 'antennas.xml')

# 5. Read coverage areas                                                    ####
cellsFileName <- file.path(path_resources, 'AntennaCells_MNO1.csv')
coverarea <- readCells(cellsFileName)
  

# 6. If method != trajectory                                                  ####
if(method != "trajectory"){
  #### ..... Read RSS                                                       ####
  
  RSS.dt <- fread(file.path(path_resources, 'RSS.csv'), 
                  colClasses = c('numeric', 'character', 'numeric', 'character', 'integer',
                                 'numeric', 'numeric', 'numeric', 'numeric', 'numeric',
                                 'integer', 'numeric'))
  
  
  
  #### ..... Set HMM prior                                                  ####
  nTiles <- gridParams$ncol * gridParams$nrow
  
  if(geolocation_prior == "uniform"){
    
    initialDistr.vec <- rep(1 / nTiles, nTiles)
    
  }
  
  if(geolocation_prior == "network"){
    
    if(emission_model == "RSS"){
      initialDistr_RSS_network.dt <- RSS.dt[
        , watt := 10**( (RSS_ori - 30) / 10 )][
        , total := sum(watt, na.rm = TRUE)][
        , list(num = sum(watt, na.rm = TRUE), total = unique(total)), by = 'rasterCell'][
        , prior_network := num / total][
        order(rasterCell)]
      initialDistr.vec <- initialDistr_RSS_network.dt$prior_network
    }
    if(emission_model == "SDM"){
      
      initialDistr_SDM_network.dt <- RSS.dt[
        , total := sum(SDM_ori, na.rm = TRUE)][
        , list(num = sum(SDM_ori, na.rm = TRUE), total = unique(total)), by = 'rasterCell'][
        , prior_network := num / total][
        order(rasterCell)]
      initialDist.vec <- initialDistr_SDM_network.dt$prior_network
    }
  }
  
  
  #### ..... Compute duplicity probabilities                                ####
  # ____ i. Get a list of detected devices                                  ####
  devices <- getDeviceIDs(events)
  
  # ____ ii. Get connections for each device                                ####
  connections <- getConnections(events)
  
  # ____ iii. Get emission probs from RSS file                              ####
  emissionProbs <- getEmissionProbs(nrows = gridParams$nrow, ncols = gridParams$ncol,
                                    signalFileName = signalFileName,
                                    sigMin = simParams$conn_threshold,
                                    handoverType = simParams$connection_type[[1]], 
                                    emissionModel = emission_model, 
                                    antennaFileName = antennasFileName)
  
  # ____ iv. Get joint emission probabilities                               ####
  jointEmissionProbs <- getEmissionProbsJointModel(emissionProbs)
  
  # ____ v. Build the generic model                                         ####
  model <- getGenericModel(
    nrows = gridParams$nrow, ncols = gridParams$ncol,
    emissionProbs, initSteady = FALSE, aprioriProb = initialDistr.vec)
  
  
  # ____ vi. Check connection-emissionProb compatibility                    ####
  check1 <- checkConnections_step1(connections, emissionProbs)
  check1$infoCheck_step1
  connections_Imp <- check1$connectionsImp
  
  # ____ vii. Check jump-emissionProb compatibility fro time_padding        ####
  check2 <- checkConnections_step2(emissionProbs = emissionProbs, 
                                   connections = connections, 
                                   gridParams = gridParams)
  connections_Imp <- check2$connections_pad
  
  
  # ____ viii. Fit models                                                   ####
  ll <- fitModels(length(devices), model, connections_Imp)
  
  # ____ ix. Build joint model                                              ####
  modelJ <- getJointModel(
    nrows = gridParams$nrow, ncols = gridParams$ncol, 
    jointEmissionProbs = jointEmissionProbs, 
    initSteady = FALSE, aprioriJointProb = initialDistr.vec)
}

# 7a. method == 1to1                                                        ####
if(method == "1to1"){
  #### ..... Set duplicity prior                                            ####
  Pii <- aprioriOneDeviceProb(simParams$prob_sec_mobile_phone, length(devices))
  
  #### ..... Build matrix of pairs of devices                               ####
  pairs4dup <- computePairs(connections_Imp, length(devices), oneToOne = TRUE)
  
  #### ..... Compute duplicity probabilities                                ####
  if(is.na(lambda)){
    
    duplicity.dt <- computeDuplicityBayesian("1to1", devices, pairs4dup, modelJ, 
                                                 ll, P1 = NULL, Pii=Pii)
  }
  if(!is.na(lambda)){
    
    duplicity.dt <- computeDuplicityBayesian("1to1", devices, pairs4dup, modelJ, 
                                                             ll, P1 = NULL, Pii=NULL,
                                                             lambda = lambda)
  }
} #end if method 1to1

# 7b. method == pairs                                                       ####
if(method == "pairs"){
  
  #### ..... Set duplicity prior                                            ####
  P1 <- aprioriDuplicityProb(simParams$prob_sec_mobile_phone, length(devices))
  
  #### ..... Build matrix of pairs of devices                               ####
  antennaNeigh <- antennaNeighbours(coverarea)
  pairs4dup<-computePairs(connections, length(devices), 
                          oneToOne = FALSE, P1=P1, limit = 0.05, 
                          antennaNeighbors = antennaNeigh)
  
  #### ..... Compute duplicity probabilities                                ####
  duplicity.dt <- computeDuplicityBayesian("pairs", devices, pairs4dup, 
                                                            modelJ, ll, P1)
  
}

# 7c. method == trajectory                                                  ####
if(method == "trajectory"){
  

  #### ..... Set duplicity prior                                            ####
  P1 <- aprioriDuplicityProb(simParams$prob_sec_mobile_phone, length(devices))
  
  #### ..... Build matrix of pairs of devices                               ####
  times <- sort(unique(events[, time]))
  antennaNeigh <- antennaNeighbours(coverarea)
  pairs4dup <- computePairs(connections, length(devices), 
                          oneToOne = FALSE, P1=P1, limit = 0.05, 
                          antennaNeighbors = antennaNeigh)
  
  #### ..... Compute duplicity probabilities                                ####
  duplicity.dt <- computeDuplicityTrajectory(path = path_postLoc,
                                             prefix = 'postLocProb_HMM_RSS_uniform',
                                             devices = devices,
                                             gridParams = gridParams,
                                             pairs = pairs4dup,
                                             P1 = P1,
                                             T = times,
                                             gamma = 0.5)
  
  
}


# 8. Treat devices with unique same connection                              ####
# all devices with the same events as another device
devIDs_sameEventVar <- devices[which(duplicated(connections) | duplicated(connections, fromLast=TRUE))]
# This will be integrated in the deduplication package
duplicity_FNcorrected.dt <- copy(duplicity.dt)[
  deviceID %in% devIDs_sameEventVar, dupP := 1]


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                              SAVE RESULTS                              ####
fwrite(duplicity.dt[, .(deviceID, dupP)],
       file.path(path_dupProb, 
                 paste0('duplicity_', method, '_', emission_model,'-', geolocation_prior, '.csv')),
       row.names = FALSE, col.names = FALSE, sep = ',')

# generar estos objetos con los eventos directamente, si misma antena dupProb=1 en ambos
# separar todo lo que usa persons a otro script de analisis

fwrite(duplicity_FNcorrected.dt[, .(deviceID, dupP)], 
       file.path(path_dupProb, 
                 paste0('duplicity_', method, '_', emission_model,'-', geolocation_prior, '-FNcorrected.csv')),
       row.names = FALSE, col.names = FALSE, sep = ',')
####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####


