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
path_metrics      <- file.path(path_root, 'metrics')
path_Nnet         <- file.path(path_root, 'data/NnetDistrib')

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)
library(destim)
library(xml2)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)              # graphics
library(latex2exp)
library(pROC)                 # for ROC curve for duplicity test
library(extraDistr)
library(bayestestR)
library(raster)               # spatial data
library(rgeos)                # spatial data
library(deduplication)
library(aggregation)

source(file.path(path_source, 'get_simConfig.R'))
source(file.path(path_source, 'get_simScenario.R'))


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                     LOAD INPUT DATA AND PARAMETERS                     ####
parameters.xml <- as_list(read_xml(file.path(path_processParam, "parameters_process.xml")))
emission_model    <- parameters.xml$process_parameters$geolocation$emission_model[[1]]
geolocation_prior <- parameters.xml$process_parameters$geolocation$prior[[1]]
method_dup        <- parameters.xml$process_parameters$deduplication$method[[1]]
nSim              <- as.integer(parameters.xml$process_parameters$aggregation$nSim[[1]])
prefix_postLoc    <- parameters.xml$process_parameters$aggregation$prefix_postLoc[[1]]

# time
simConfigParam.list <- get_simConfig(path_simConfig)
initialTime <- as.integer(simConfigParam.list$simulation.xml$simulation$start_time)
finalTime   <- as.integer(simConfigParam.list$simulation.xml$simulation$end_time)
t_increment <- as.integer(simConfigParam.list$simulation.xml$simulation$time_increment)
times       <- seq(from = initialTime, to = (finalTime - t_increment), by = t_increment)

# File names
dpFileName <-   file.path(
  path_dupProb, 
  paste0('duplicity_', method_dup, '_', emission_model,'-', geolocation_prior, '-FNcorrected.csv'))

NnetFileName <- file.path(path_Nnet, paste0('Nnet_', emission_model, '-', geolocation_prior, '_', method_dup, '.csv'))

fileGridName <- file.path(path_resources, simConfigParam.list$filesNames$fileGridName)

rgFileName <- file.path(path_resources, 'regions.csv')






####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####          ESTIMATE NUMBER OF INDIVIDUALS DETECTED BY THE NETWORK        ####
Nnet.dt <- rNnetEvent(n = nSim, 
                      gridFileName = fileGridName, 
                      dupFileName = dpFileName, 
                      regsFileName = rgFileName, 
                      postLocPath = path_postLoc, 
                      prefix = prefix_postLoc, 
                      times = times)


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####          SAVE NUMBER OF INDIVIDUALS DETECTED BY THE NETWORK            ####
fwrite(Nnet.dt, NnetFileName)

