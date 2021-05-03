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
library(tmap)
library(mobvis)

# Function get_simConfig to read the input files of the simulator
source(file.path(path_source, 'get_simConfig.R'))
# Function get_simScenario to read the output files of the simulator
source(file.path(path_source, 'get_simScenario.R'))
# Function tileEquivalence to compute the equivalence between rastercell (R) and tiles (simulator)
source(file.path(path_source, 'tileEquivalence.R'))
# Function to fit and compute HMM model with the events of a specific device
source(file.path(path_source, 'compute_HMM.R'))
# Function to transform de output of compute_HMM
source(file.path(path_source, 'transform_postLoc.R'))
# Function to fit static model with uniform and network priors
source(file.path(path_source, 'compute_staticModel.R'))



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####                         LOAD PARAMETERS                              #####
sim <- list(simConfig_dir = path_simConfig,
            resources_dir = path_resources,
            networkEvents = path_events,
            mno = "MNO1",
            crs = sf::st_crs(NA))

# cell to draw
cellID <- 29


rst <- sim_get_raster(sim)
cp  <- sim_get_cellplan(sim)
map <- sim_get_region(sim)
strength <- sim_get_signal_strength(sim, rst, cp)


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####                   PLOT RSS, SDM, BEST SERVER                         #####
map_sig_strength(rst,
                 dt = strength,
                 cp = cp,
                 region = map,
                 cells = cellID,
                 type = "dBm",
                 interactive = FALSE,
                 settings = mobvis_settings(cell_size = 1,
                                            cell_labels = TRUE))

map_sig_strength(rst,
                 dt = strength,
                 cp = cp,
                 region = map,
                 cells = cellID,
                 type = "s",
                 interactive = FALSE,
                 settings = mobvis_settings(cell_size = 1,
                                            cell_labels = TRUE))

map_best_server(rst = rst,
                dt = strength,
                cp = cp,
                region = region,
                cells = cellID,
                type = "bsm",
                interactive = FALSE,
                settings = mobvis_settings(cell_size = 1,
                                           cell_labels = TRUE))











