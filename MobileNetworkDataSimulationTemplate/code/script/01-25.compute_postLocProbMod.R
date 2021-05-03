####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
path_root         <- 'D:/Documentos/Universidad/TFG/MobileNetworkDataSimulationTemplate'
path_source       <- file.path(path_root, 'code/src')
path_simConfig    <- file.path(path_root, 'data/simulatorConfig')
path_events       <- file.path(path_root, 'data/networkEvents')
path_eventLoc     <- file.path(path_root, 'data/eventLocProb')
path_resources    <- file.path(path_root, 'param/resources')
path_processParam <- file.path(path_root, 'param/process')
path_postLoc      <- file.path(path_root, 'data/postLocProb')



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)           # manage data
library(xml2)                 # to read simulatorConfig and parameters
library(tidyr)                # transform xml documents to tables (tibbles)
library(rgeos)                # spatial data
library(destim)               # HMM model
library(stringr)              # to pad strings
library(Matrix)
library(Rsolnp)
library(FactoMineR)
library(factoextra)
library(ISLR)
library(flashClust)
library(gridExtra)
library(dplyr)

# Function get_simConfig to read the input files of the simulator
source(file.path(path_source, 'get_simConfig.R'))
# Function get_simScenario to read the output files of the simulator
source(file.path(path_source, 'get_simScenario.R'))
# Function tileEquivalence to compute the equivalence between rastercell (R) and tiles (simulator)
source(file.path(path_source, 'tileEquivalence.R'))
# Function to fit and compute HMM model with the events of a specific device
source(file.path(path_source, 'compute_HMMParams.R'))
source(file.path(path_source, 'compute_HMMPostLoc.R'))
source(file.path(path_source, 'compute_HMM.R'))
# Function to transform de output of compute_HMM
source(file.path(path_source, 'transform_postLoc.R'))
# Function to fit static model with uniform and network priors
source(file.path(path_source, 'compute_staticModel.R'))


#### :::::::::::::::::::::::::::::::::::::: ####
#####       START OF EXPERIMENTS           #####


#Cambiar longitud del bucle al gusto para hacer los experimentos deseados, lo hace sobre las columnas de
#parameters.method, csv creado manualmente con los parametros de cada experimento

parameters.method <- fread(file.path(path_processParam, "parameters_method.csv"))

for(m in 1:dim(parameters.method)[1]){
  
  
####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####           LOAD PARAMETERS AND NETWORK EVENT DATA                     #####

parameters.xml    <- as_list(read_xml(file.path(path_processParam, "parameters_process.xml")))
geolocation_model <- parameters.xml$process_parameters$geolocation$model_name[[1]]
geolocation_prior <- parameters.xml$process_parameters$geolocation$prior[[1]]
emission_model    <- parameters.xml$process_parameters$geolocation$emission_model[[1]]
transition_model    <- parameters.xml$process_parameters$geolocation$transition_model[[1]]

simConfigParam.list <- get_simConfig(path_simConfig)

# time
times       <- simConfigParam.list$simConfigParameters$times
t_increment <- simConfigParam.list$simConfigParameters$t_increment

# grid
fileGridName       <- file.path(path_resources, simConfigParam.list$filesNames$fileGridName)
fileEventsInfoName <- file.path(path_events, simConfigParam.list$filesNames$fileEventsInfoName)
simScenario.list   <- get_simScenario(fileGridName = fileGridName, fileEventsInfoName = fileEventsInfoName)

gridParam  <- simScenario.list$gridParam
ncol_grid  <- gridParam[['No Tiles Y']]
nrow_grid  <- gridParam[['No Tiles X']]
tile_sizeX <- gridParam[['X Tile Dim']]
tile_sizeY <- gridParam[['Y Tile Dim']]
ntiles_x   <- gridParam[['No Tiles X']]
ntiles_y   <- gridParam[['No Tiles Y']]
ntiles     <- ntiles_x * ntiles_y

# tile-rasterCell equivalence
tileEquiv.dt <- simScenario.list$tileEquiv.dt

# RSS and other network parameters
RSS.dt <- fread(file.path(path_resources, 'RSS.csv'), 
                colClasses = c('numeric', 'character', 'numeric', 'character', 'integer',
                               'numeric', 'numeric', 'numeric', 'numeric', 'numeric',
                               'integer', 'numeric'))

# Maximum velocity
vMax_ms <- as.numeric(parameters.xml$process_parameters$geolocation$params$vmax_ms[[1]])

# Network event data
events.dt    <- simScenario.list$events.dt


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                      SET TIME PADDING PARAMETERS                       ####
distMax <- vMax_ms * t_increment
pad_coef <- as.integer(ceiling(distMax / max(tile_sizeX, tile_sizeY)))
pad_coef <- pad_coef + 1 



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                INITIAL STATE DISTRIBUTION (PRIOR)                      ####
# Prepare prior_network distribution

if(geolocation_prior == "uniform"){
  
  prior_network <- rep(1 / ntiles, ntiles)
  
}

if(geolocation_prior == "network"){
  
  if(emission_model == "RSS"){
    
    initialDistr_RSS_network.dt <- RSS.dt[
      , watt := 10**( (RSS_ori - 30) / 10 )][
      , total := sum(watt, na.rm = TRUE)][
      , list(num = sum(watt, na.rm = TRUE), total = unique(total)), by = 'rasterCell'][
      , prior_network := num / total][
      order(rasterCell)]
    prior_network <- initialDistr_RSS_network.dt$prior_network
    
  }
  
  if(emission_model == "SDM"){
    
    initialDistr_SDM_network.dt <- RSS.dt[
      , total := sum(SDM_ori, na.rm = TRUE)][
      , list(num = sum(SDM_ori, na.rm = TRUE), total = unique(total)), by = 'rasterCell'][
      , prior_network := num / total][
      order(rasterCell)]
    prior_network <- initialDistr_SDM_network.dt$prior_network
    
  }
}

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                          EMISSION MODEL                                ####
# Load the event location probabilities:
# (i)  For the HMM this matrix contains the emission probabilities
# (ii) For the static analyses this matrix contains the likelihoods to apply Bayes' rule
emissionProb_rasterCell.dt <- fread(file.path(path_eventLoc, "eventLocProb.csv"),
                                    colClasses = c('character', 'numeric', 'numeric', 'character', 'numeric'),
                                    sep = ',')[
  tileEquiv.dt, on = 'tile'][
  , c('device', 'time', 'tile') := NULL][
  order(event_cellID, rasterCell)]
setcolorder(emissionProb_rasterCell.dt, c('rasterCell', 'event_cellID', 'eventLocProb'))
emissionProb_rasterCell.matrix <- as.matrix(
  dcast(emissionProb_rasterCell.dt, 
        rasterCell ~ event_cellID, value.var = 'eventLocProb')[
  , rasterCell := NULL])


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####                          TRANSITION MODEL                            #####

if(transition_model == "HMMrectangle"){
  
  ## Initialize HMM                                                           ####
  cat('Initializing HMM...')
  model <- HMMrectangle(nrow_grid, ncol_grid)
  emissions(model) <- emissionProb_rasterCell.matrix # eventLoc for each antenna
  
  model <- initparams(model)  # initialize transition prob
  model <- minparams(model)   # parameter reduction according to restrictions
  
  istates(model) <- prior_network 
  
  cat('ok.\n\n')
  
}


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####     COMPUTE POSTERIOR LOCATION PROBABILITIES AND INDICATORS           #####
cat('Computing posterior location probabilities...\n')

postLocProb.list      <- list()
postLocJointProb.list <- list()

deviceIDs <- sort(unique(events.dt$device))
RSS.dt <- merge(RSS.dt, emissionProb_rasterCell.dt[, .(rasterCell, event_cellID, eventLocProb)], 
                by.x = c('rasterCell', 'antennaID'), by.y = c('rasterCell', 'event_cellID'),
                all.x = TRUE)




#### :::::::::::::::::::::::::::::::::::::: ####
#####       START OF EXPERIMENTS           #####


#Cambiar longitud del bucle al gusto para hacer los experimentos deseados, lo hace sobre las columnas de
#parameters.method, csv creado manualmente con los parametros de cada experimento


  #### :::::::::::::::::::::::::::::::::::::: ####
  #####     COMPUTING KEY PARAMETER          #####
 
  
  
  if (parameters.method[m,1] == 'mean'){
    
    cat(paste0('Computing Key parameter with case ', parameters.method[m,2],'...\n'))  
    
    n=as.integer(as.numeric(parameters.method[m,3])*length(deviceIDs))
    
    deviceSample <- sample(deviceIDs, n, replace=F)
    ParameterMatrix <- matrix(nrow = n, ncol = 2)
  
  
    
    for(i in seq(along = deviceSample)){
      
      ###                  :: For each device                          ####
      devID <- deviceSample[i]
      cat(paste0('    device ', devID,'...\n'))
      
      # Selecting device events
      cat('       selecting network events...')
      events_device.dt <- events.dt[
        device == devID, .(device, time, antennaID)][
          order(device, time)]
      
      antennas_deviceID  <- unlist(events_device.dt[, c("antennaID")])
      cat(' ok.\n')  
      
      if(!all(is.na(antennas_deviceID))){
        HMMmodel_outputParam <- compute_HMMParams(model = model, observedValues = antennas_deviceID, 
                                                  pad_coef = pad_coef, init = TRUE)
        ParameterMatrix[i,] <- HMMmodel_outputParam$parameters$reducedparams$params
        
      }
    }
    
    FinalParameters=c(mean(ParameterMatrix[,1]), mean(ParameterMatrix[,2]))
    model$parameters$reducedparams$params = FinalParameters
    model$parameters$transitions=gettransmatrix(model) %*% c(rparams(model), 1) 
  
  }
  
  
  if (parameters.method[m,1] == 'cluster'){
    cat(paste0('Computing Key parameter with case ', parameters.method[m,2],'...\n'))  
    antennascoord <- matrix(nrow = 70, ncol = 2)
    masterdevices <- matrix(nrow = length(deviceIDs), ncol = 2*(length(simConfigParam.list[["simConfigParameters"]][["times"]])))
    antennascoord[,1] <- simConfigParam.list[["antennas_xml"]][["x"]]
    antennascoord[,2] <- simConfigParam.list[["antennas_xml"]][["y"]]
    
    
    for(i in seq(along = deviceIDs)){
      
      ###                  :: For each device                          ####
      devID <- deviceIDs[i]
      #cat(paste0('    device ', devID,'...\n'))
      
      # Selecting device events
      #cat('       selecting network events...')
      events_device.dt <- events.dt[
        device == devID, .(device, time, antennaID)][
          order(device, time)]
      
      antennas_deviceID  <- unlist(events_device.dt[, c("antennaID")])
    #cat(' ok.\n')  
      
      antennas_deviceIDcoord=NA
      
      for(f in seq(along = antennas_deviceID)){
        antennas_deviceIDcoord=c(antennas_deviceIDcoord,antennascoord[as.numeric(antennas_deviceID[f]),1],antennascoord[as.numeric(antennas_deviceID[f]),2])
      }
      
      antennas_deviceIDcoord=antennas_deviceIDcoord[-1]
      
      masterdevices[i,]= antennas_deviceIDcoord
    }
    
    col_index_PAR = seq(2,ncol(masterdevices), by=2)
    col_index_IMPAR = seq(1,ncol(masterdevices), by=2)
    
    masterdevicesX = masterdevices[,col_index_IMPAR]
    masterdevicesY = masterdevices[,col_index_PAR]
     
      
    finalVectorsX <- t(apply(masterdevicesX,1,diff))
    finalVectorsY <- t(apply(masterdevicesY,1,diff))
    
    finalvectors=matrix(nrow=nrow(masterdevices), ncol=2*ncol(finalVectorsX))
    
    col_index_PAR = seq(2,ncol(finalvectors), by=2)
    col_index_IMPAR = seq(1,ncol(finalvectors), by=2)
    
    finalvectors[,col_index_IMPAR]= finalVectorsX
    finalvectors[,col_index_PAR]= finalVectorsY
    
    dist_eucl <- dist(finalvectors, method = 'euclidean')
    
    
    #Hice kmeans sobre masterdevices para comprobar que da mejores resultados que sobre finalvectors
    
    Ncluster=as.integer(parameters.method[m,3])
    kmeans_cluster <- kmeans( x = dist_eucl, centers = Ncluster)
    # kmeans_cluster[["centers"]]
    fviz_nbclust(as.matrix(dist_eucl), kmeans, method = "wss") +
           geom_vline(xintercept = 5, linetype = 2) +
           theme_bw()
    
    clusterdevice <- kmeans_cluster[["cluster"]]
    devicecluster <- matrix(nrow = Ncluster, ncol = length(deviceIDs))
    for (b in seq(along = clusterdevice)){
      p=1
      while (!is.na(devicecluster[clusterdevice[b],p])){
        p=p+1
      }
      devicecluster[clusterdevice[b],p]=b
    }
    
    ParameterMatrix <- matrix(nrow = Ncluster, ncol = 2) 
    
    for (i in seq(Ncluster)){
      number_of_nas <- sum(is.na(devicecluster[i,]))
      devposition <- sample(x = devicecluster[i,-((length(deviceIDs)-number_of_nas+1):length(deviceIDs))], size = 1)
      
      devID <- deviceIDs[devposition]
      #(paste0('    device ', devID,'...\n'))
      
      # Selecting device events
      #cat('       selecting network events...')
      events_device.dt <- events.dt[
        device == devID, .(device, time, antennaID)][
          order(device, time)]
      
      antennas_deviceID  <- unlist(events_device.dt[, c("antennaID")])
      #cat(' ok.\n')  
      
      HMMmodel_outputParam <- compute_HMMParams(model = model, observedValues = antennas_deviceID, 
                                                pad_coef = pad_coef, init = TRUE)
      ParameterMatrix[i,] <- HMMmodel_outputParam$parameters$reducedparams$params
      
    }
      
  }
    
    
  
  
  
  ####  :::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####    COMPUTING POSTERIOR LOCATION PROBABILITIES     #####
  
  for(i in seq(along = deviceIDs)){
    
    ###                  :: For each device                          ####
    devID <- deviceIDs[i]
    cat(paste0('    device ', devID,'...\n'))
    
    # Selecting device events
    cat('       selecting network events...')
    events_device.dt <- events.dt[
      device == devID, .(device, time, antennaID)][
      order(device, time)]
  
    antennas_deviceID  <- unlist(events_device.dt[, c("antennaID")])
    cat(' ok.\n')  
    
    if(!all(is.na(antennas_deviceID))){
  
      #####                Fit and compute HMM model                 ####
      
      if (parameters.method[m,1] == 'cluster'){
        cat(paste0('    cluster ', clusterdevice[i],'...\n'))
        model$parameters$reducedparams$params = ParameterMatrix[clusterdevice[i],]
        model$parameters$transitions=gettransmatrix(model) %*% c(rparams(model), 1) 
      }
      
      if(parameters.method[m,1] == 'All'){
        cat('       fit and compute HMM model...\n')
        HMMmodel_output <- compute_HMM(model = model, observedValues = antennas_deviceID, 
                                       pad_coef = pad_coef, init = TRUE)
          
      }else{
      
        cat('compute HMM model...\n')
        HMMmodel_output <- compute_HMMPostLoc(model = model, observedValues = antennas_deviceID, 
                                     pad_coef = pad_coef)
      }
      # model_devID <- HMMmodel_output$model_devID
      postLocProb_HMM_deviceID.matrix      <- HMMmodel_output$postLocP 
      postLocJointProb_HMM_deviceID.matrix <- HMMmodel_output$postJointLocP
      
      cat(' ok.\n')
      #####                 Transform output HMM model                 ####
      cat('       transform output HMM model...\n')
      transform_output <- transform_postLoc(postLocP = postLocProb_HMM_deviceID.matrix,
                                           postLocJointP = postLocJointProb_HMM_deviceID.matrix,
                                           observedValues = antennas_deviceID,
                                           times = times, t_increment = t_increment, 
                                           ntiles = ntiles, pad_coef = pad_coef,
                                           tileEquiv.dt = tileEquiv.dt, devID = devID,
                                           sparse_postLocP = TRUE, sparse_postLocJointP = TRUE)
      rm(postLocProb_HMM_deviceID.matrix)
      rm(postLocJointProb_HMM_deviceID.matrix)
      gc()
      postLocProb_deviceID.dt <- transform_output$postLocProb
      postLocProb.list[[as.character(devID)]] <- postLocProb_deviceID.dt[
        , c('device', 'time', 'tile', 'event_cellID', 'postLocProb'), with = FALSE]
      
      postLocJointProb_deviceID.dt <- transform_output$postLocJointProb
      postLocJointProb.list[[as.character(devID)]] <- postLocJointProb_deviceID.dt[
        , time_to := time_from + t_increment][
        , c('device', 'time_from', 'time_to', 'tile_from', 'tile_to', 'event_cellID_from', 'event_cellID_to', 'postLocProb'), with = FALSE]
      
      # fwrite(postLocProb_deviceID.dt[, .(device, time, tile, event_cellID, postLocProb)], 
      #        file.path(path_postLoc, 
      #                  paste0('postLocProb_', geolocation_model, '_', emission_model,
      #                         '-', geolocation_prior, '_', devID, '.csv')),
      #        col.names = FALSE, row.names = FALSE, sep = ',')
      # 
      # fwrite(postLocJointProb_deviceID.dt[, .(device, time_from, time_to, tile_from, tile_to, event_cellID_from, event_cellID_to, postLocProb)], 
      #        file.path(path_postLoc, 
      #                  paste0('postLocJointProb_', geolocation_model, '_', emission_model,
      #                         '-', geolocation_prior, '_', devID, '.csv')),
      #         col.names = FALSE, row.names = FALSE, sep = ',')
      
      rm(postLocJointProb_deviceID.dt)
      rm(transform_output)
      gc()
      cat(' ok.\n')  
  
      
    }
    ###                  :: End For each device                          ####
    
  }  
  
  
  postLocProb.dt      <- rbindlist(postLocProb.list)
  if(geolocation_model == "HMM"){
    postLocJointProb.dt <- rbindlist(postLocJointProb.list)
  }
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                   SAVE POSTERIOR LOCATION PROBS                      #####
  fwrite(postLocProb.dt, 
         file.path(path_postLoc, 
                   paste0('postLocProb_', geolocation_model, '_', 
                          emission_model, '-', geolocation_prior, '_', parameters.method[m,2], '.csv')))
  
  if(geolocation_model == "HMM"){
    
    fwrite(postLocJointProb.dt, 
           file.path(path_postLoc, 
                     paste0('postLocJointProb_', geolocation_model, '_', 
                            emission_model, '-', geolocation_prior, '_', parameters.method[m,2], '.csv')))
  }

  
}