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
path_Ntarget      <- file.path(path_root, 'data/NtargetDistrib')

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)
library(destim)
library(stringr)
library(xml2)
library(tibble)
library(tidyr)
library(ggplot2)              # graphics
library(latex2exp)
library(rgeos)
library(pROC)                 # for ROC curve for duplicity test
library(extraDistr)
library(bayestestR)
library(raster)               # spatial data
library(deduplication)
library(aggregation)
library(inference)

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
inference_model   <- parameters.xml$process_parameters$inference$model_name[[1]]

dpFileName <-   file.path(
  path_dupProb, 
  paste0('duplicity_', method_dup, '_', emission_model,'-', geolocation_prior, '-FNcorrected.csv'))

NnetFileName       <- file.path(path_Nnet, paste0('Nnet_', emission_model, '-', geolocation_prior, '_', method_dup, '.csv'))
Nnet.dt            <- readNnetInitial(NnetFileName)
NnetODFileName     <- file.path(path_Nnet, paste0('NnetOD_', emission_model, '-', geolocation_prior, '_', method_dup, '.csv'))
NnetODFileName_zip <- gsub('.csv', '.zip', NnetODFileName)


rgFileName        <- file.path(path_resources, 'regions.csv')
popRegFileName    <- file.path(path_resources, 'popRegister.csv')
pentrRateFileName <- file.path(path_resources, 'pentrRate.csv')
gridFileName      <- file.path(path_resources, 'grid.csv')

NtargetFileName        <- file.path(path_Ntarget, paste0('Ntarget_', emission_model, '-', geolocation_prior, '_', method_dup, '_', inference_model, '.csv'))
NtargetODFileName      <- file.path(path_Ntarget, paste0('NtargetOD_', emission_model, '-', geolocation_prior, '_', method_dup, '_', inference_model, '.csv'))
statsNtargetFileName   <- file.path(path_Ntarget, paste0('statsNtarget_', emission_model, '-', geolocation_prior, '_', method_dup, '_', inference_model, '.csv'))
statsNtargetODFileName <- file.path(path_Ntarget, paste0('statsNtargetOD_', emission_model, '-', geolocation_prior, '_', method_dup, '_', inference_model, '.csv'))


# time
simConfigParam.list <- get_simConfig(path_simConfig)
initialTime <- as.integer(simConfigParam.list$simulation.xml$simulation$start_time)
finalTime   <- as.integer(simConfigParam.list$simulation.xml$simulation$end_time)
t_increment <- as.integer(simConfigParam.list$simulation.xml$simulation$time_increment)
times       <- seq(from = initialTime, to = (finalTime - t_increment), by = t_increment)



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                     COMPUTE DEDUPLICATION FACTORS                      ####
omega_r <- computeDeduplicationFactors(dpFileName, rgFileName, prefix_postLoc, path_postLoc)

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                     COMPUTE MODEL (HYPER)PARAMETERS                    ####
params <- computeDistrParams(omega_r, popRegFileName, pentrRateFileName, rgFileName, gridFileName)

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                        COMPUTE INITIAL POPULATION                      ####
N_initial <- computeInitialPopulation(Nnet.dt, params, popDistr = inference_model, rndVal = TRUE)
  

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                        COMPUTE POPULATION AT t > t0                    ####
N_t <- computePopulationT(N_initial$rnd_values, NnetODFileName_zip, rndVal = TRUE)

colNames <- names(N_t[[1]][['stats']])
stats_N_t <- lapply(N_t, `[[`, 'stats')
stats_N_t <- rbindlist(lapply(times, function(t){
  out <- stats_N_t[[as.character(t)]][, time:= t]}))
setcolorder(stats_N_t, c('time', colNames))

colNames <- names(N_t[[1]][['rnd_values']])
N_t <- lapply(N_t, `[[`, 'rnd_values')
N_t <- rbindlist(lapply(times, function(t){
  out <- N_t[[as.character(t)]][, time:= t]}))
setcolorder(N_t, unique(c('time', colNames)))

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           COMPUTE OD MATRICES                          ####
OD_t <- computePopulationOD(N_initial$rnd_values, NnetODFileName_zip, rndVal = TRUE)

colNames <- names(OD_t[[1]][['stats']])
stats_OD_t <- rbindlist(mapply(set, x = lapply(OD_t_cp, `[[`, 'stats'), j = 'time_from', value = times[-length(times)], SIMPLIFY = FALSE))[
  , time_to := time_from + t_increment]
setcolorder(stats_OD_t, c('time_from', 'time_to', colNames))


colNames <- names(OD_t[[1]][['rnd_values']])
OD_t <- rbindlist(mapply(set, x = lapply(OD_t, `[[`, 'rnd_values'), j = 'time_from', value = times[-length(times)], SIMPLIFY = FALSE))[
  , time_to := time_from + t_increment]
setcolorder(OD_t, c('time_from', 'time_to', colNames))


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                              SAVE OUTPUTS                              ####
fwrite(N_t, NtargetFileName)
fwrite(OD_t, NtargetODFileName)
fwrite(stats_N_t, statsNtargetFileName)
fwrite(stats_OD_t, statsNtargetODFileName)