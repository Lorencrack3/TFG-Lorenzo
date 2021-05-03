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
ci_prob           <- as.numeric(parameters.xml$process_parameters$inference$ci_prob[[1]])
ci_method         <- parameters.xml$process_parameters$inference$ci_method[[1]]


# time
simConfigParam.list <- get_simConfig(path_simConfig)
initialTime <- as.integer(simConfigParam.list$simulation.xml$simulation$start_time)
finalTime   <- as.integer(simConfigParam.list$simulation.xml$simulation$end_time)
t_increment <- as.integer(simConfigParam.list$simulation.xml$simulation$time_increment)
times       <- seq(from = initialTime, to = (finalTime - t_increment), by = t_increment)


# Filenames
dpFileName <-   file.path(
  path_dupProb, 
  paste0('duplicity_', method_dup, '_', emission_model,'-', geolocation_prior, '-FNcorrected.csv'))

NnetFileName       <- file.path(path_Nnet, paste0('Nnet_', emission_model, '-', geolocation_prior, '_', method_dup, '.csv'))
Nnet_initial.dt    <- readNnetInitial(NnetFileName)
setnames(Nnet_initial.dt, 'N', 'Nnet')
NnetODFileName     <- file.path(path_Nnet, paste0('NnetOD_', emission_model, '-', geolocation_prior, '_', method_dup, '.csv'))
NnetOD.dt          <- fread(NnetODFileName)


rgFileName        <- file.path(path_resources, 'regions.csv')
popRegFileName    <- file.path(path_resources, 'popRegister.csv')
pentrRateFileName <- file.path(path_resources, 'pentrRate.csv')
gridFileName      <- file.path(path_resources, 'grid.csv')

NtargetFullFileName      <- file.path(path_Ntarget, paste0('Ntarget_full_', emission_model, '-', geolocation_prior, '_', method_dup, '.csv'))
statsNtargetFullFileName <- file.path(path_Ntarget, paste0('statsNtarget_full_', emission_model, '-', geolocation_prior, '_', method_dup, '.csv'))

NtargetFullFileName        <- file.path(path_Ntarget, paste0('Ntarget_full_', emission_model, '-', geolocation_prior, '_', method_dup, '_', inference_model, '.csv'))
NtargetODFullFileName      <- file.path(path_Ntarget, paste0('NtargetOD_full_', emission_model, '-', geolocation_prior, '_', method_dup, '_', inference_model, '.csv'))
statsNtargetFullFileName   <- file.path(path_Ntarget, paste0('statsNtarget_full_', emission_model, '-', geolocation_prior, '_', method_dup, '_', inference_model, '.csv'))
statsNtargetODFullFileName <- file.path(path_Ntarget, paste0('statsNtargetOD_full_', emission_model, '-', geolocation_prior, '_', method_dup, '_', inference_model, '.csv'))

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                     COMPUTE DEDUPLICATION FACTORS                      ####
omega_r <- computeDeduplicationFactors(dpFileName, rgFileName, prefix_postLoc, path_postLoc)

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                     COMPUTE MODEL (HYPER)PARAMETERS                    ####
params <- computeDistrParams(omega_r, popRegFileName, pentrRateFileName, rgFileName, gridFileName)


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                COMPUTE POSTERIOR DISTRIBUTION                          ####

### 1. Initial time                                                         ####
if (inference_model == 'BetaNegBin') {
  
  N_full_initial <- Nnet_initial.dt[
    , row := .I][
    params, on = 'region'][
    , list(region = region,
           Nnet = Nnet,
           NPop = Nnet + rbnbinom(1, Nnet + 1, alpha - 1, beta)), by = 'row']
  
  stats_N_full_initial <- computeStats(N_full_initial, ci_prob, ci_method, by = 'region')
}


if (inference_model == 'NegBin') {
  
  N_full_initial <- Nnet_initial.dt[
    , row := .I][
      params, on = 'region'][
        , list(region = region,
               Nnet = Nnet,
               NPop = Nnet + rnbinom(1, Nnet + 1, (alpha - 1) / (alpha + beta - 1))), by = 'row']
  
  stats_N_full_initial <- computeStats(N_full_initial, ci_prob, ci_method, by = 'region')
}


if (inference_model == 'StateNegBin') {
  
  N_full_initial <- Nnet_initial.dt[
    , row := .I][
      params, on = 'region'][
        , list(region = region,
               Nnet = Nnet,
               NPop = Nnet + rnbinom(1, Nnet + zeta + 1, 1 - beta * Q / (alpha + beta))), by = 'row']
                      
  stats_N_full_initial <- computeStats(N_full_initial, ci_prob, ci_method, by = 'region')
}


#### 2. Times t > t0                                                        ####
NnetOD_totalFrom.dt <- NnetOD.dt[
  , list(Nnet_totalFrom = sum(Nnet)), by = c('time_from', 'time_to', 'region_from', 'iter')]
tauOD.dt <- merge(NnetOD.dt, NnetOD_totalFrom.dt, 
                  by = c('time_from', 'time_to', 'region_from', 'iter'))[
  , tau_Nnet := Nnet / Nnet_totalFrom][
  Nnet_totalFrom == 0, tau_Nnet := 0][
  , region_from := factor(region_from)][
  , region_to := factor(region_to)][
  , Nnet_totalFrom := NULL][
  , Nnet := NULL]

N_full_initial[
  , iter := rep(1:nSim, times = nrow(N_full_initial) / nSim)][
  , row := NULL]

N_full_t.list <- vector(mode = 'list', length = length(times) - 1)
names(N_full_t.list) <- as.character(times[-length(times)])
runningN.dt <- copy(N_full_initial)[
  , region := factor(region)][
  , Nnet := NULL]
setcolorder(runningN.dt, c('region',  'iter', 'NPop'))
N_full_t.list[[as.character(initialTime)]] <- runningN.dt[
  , time := initialTime]
for (timeFrom in times[-length(times)]){
  
  cat(paste0('Time ', timeFrom, '... '))  
  tempDT <- merge(tauOD.dt[time_from == timeFrom], runningN.dt,
                  by.x = c('region_from', 'iter'), by.y = c('region', 'iter'), 
                  allow.cartesian = TRUE)
  runningN.dt <- tempDT[, list(NPop = sum(tau_Nnet * NPop)), by = c('region_to', 'iter')][
    , time := timeFrom + t_increment]
  rm(tempDT)
  setnames(runningN.dt, c('region_to'), c('region'))
  N_full_t.list[[as.character(timeFrom + t_increment)]] <- runningN.dt
  gc()
  cat(paste0(' ok.\n'))
}

N_full_t <- rbindlist(N_full_t.list)

stats_N_full_t <- rbindlist(lapply(times, function(t){
  
  tempDT <- N_full_t.list[[as.character(t)]]
  out <- computeStats(tempDT, ciprob = ci_prob, method = ci_method, by = 'region')[
    , time := as.integer(t)]
  return(out)
}))

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                              SAVE OUTPUTS                              ####
fwrite(N_full_t, NtargetFullFileName)
fwrite(stats_N_full_t, statsNtargetFullFileName)
