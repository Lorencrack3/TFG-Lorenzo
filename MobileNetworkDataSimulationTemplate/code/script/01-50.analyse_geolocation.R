####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
path_root         <- 'D:/Github/TFG-Lorenzo/MobileNetworkDataSimulationTemplate'
path_source       <- file.path(path_root, 'code/src')
path_simConfig    <- file.path(path_root, 'data/simulatorConfig')
path_events       <- file.path(path_root, 'data/networkEvents')
path_eventLoc     <- file.path(path_root, 'data/eventLocProb')
path_resources    <- file.path(path_root, 'param/resources')
path_processParam <- file.path(path_root, 'param/process')
path_postLoc      <- file.path(path_root, 'data/postLocProb')
path_experiments  <- file.path(path_root, 'data/postLocProbExperiments')
path_img          <- file.path(path_root, 'metrics/img')
path_grTruth      <- file.path(path_root, 'param/groundTruth')


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)           # manage data
library(tmap)
library(mobvis)
library(mobloc)
library(xml2)
library(magrittr)
library(tibble)
library(tidyr)
library(stringr)
library(rgeos)
library(ggplot2)
library(latex2exp)
library(MASS)
library(ISLR)
library(dplyr)

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
# Function to compute center of probabilities (to be merged with that from deduplication)
source(file.path(path_source, 'cp.R'))
# Function to compute radius of dispersion
source(file.path(path_source, 'rd.R'))


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                      READ RDS DATA FILES                               ####
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
tile_size  <- mean(c(tile_sizeX, tile_sizeY))
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

# Network event data
events.dt    <- simScenario.list$events.dt

# Centroid coordinates of each tile
centroidCoord.dt <- fread(file.path(path_resources, 'centroidCoord.csv'))

# True positions of each device
truePositions_device.dt <- fread(file.path(path_grTruth, 'truePositions_device.csv'))


msdmaster.dt=list()
parameters.method <- fread(file.path(path_processParam, "parameters_method.csv"))



#En este bucle calculamos el sesgo, rmsd y msd de todos los experimentos y guardamos los 
#resultados en variables internas

for(m in 1:dim(parameters.method)[1]){
# for(m in 1:1){
  postLocProbAux.dt <- fread(file.path(path_postLoc, paste0('postLocProb_', geolocation_model, '_RSS-', geolocation_prior, '_', parameters.method[m,2], '.csv')))
  postLocProbAux.dt <- merge(postLocProbAux.dt, centroidCoord.dt, by = 'tile') 
  
  #### ** Bias                                                                ####
  cp_truePosAux.dt <- postLocProbAux.dt[
    , as.list(cp(x = centroidCoord_x, y = centroidCoord_y, w = postLocProb)), by = c('device', 'time')][
      truePositions_device.dt[, .(time, device, x, y)], on = c('device', 'time')][
        , cp_tp := sqrt( (cp_x - x)**2 + (cp_y - y)**2)][
          , .(device, time, cp_tp)][
            order(device, time)]
  
  #### ** rmsd                                                                ####
  rmsdAux.dt <- postLocProbAux.dt[
    , list(rd = rd(x = centroidCoord_x, y = centroidCoord_y, w = postLocProb)), by = c('device', 'time', 'cluster', 'event_cellID')]
  
  msdAux.dt <- merge(cp_truePosAux.dt, rmsdAux.dt, by = c('device', 'time'))

  
  msdAux.dt[
    , msd := cp_tp**2 + rd**2][
      , model := paste0(geolocation_model, '-', geolocation_prior)][
        , experiment := parameters.method[m,2]
      ]
  
  #Añadimos el numero de cambios de antena por dispositivo a msdAux.
  
  antenna_changes = msdAux.dt %>% 
    group_by(device) %>%
    summarise(antenna_changes = length(unique(event_cellID)))
  
  msdAux.dt=merge(msdAux.dt, antenna_changes, by = 'device')
  
  #MSD maestro con todos los experimentos
  msdmaster.dt <- rbindlist(list(msdmaster.dt,msdAux.dt), fill = TRUE)
  #msdmaster.dt <- msdmaster.dt[,as.list(cp(x = centroidCoord_x, y = centroidCoord_y, w = postLocProb)), by = c('device', 'time', 'cluster')]
  
  
  #Asignaciones a variables internas para facilitar su posterior analisis
  
  masterpostloc=paste0('postLocProb_', parameters.method[m,2])
  cp_truePos=paste0('cp_truePos_', parameters.method[m,2])
  rmsd=paste0('rmsd_', parameters.method[m,2])
  msd=paste0('msd_', parameters.method[m,2])
  assign(masterpostloc, fread(file.path(path_postLoc, paste0('postLocProb_', geolocation_model, '_RSS-', geolocation_prior, '_', parameters.method[m,2], '.csv'))) )
  assign(cp_truePos, cp_truePosAux.dt)
  assign(rmsd, rmsdAux.dt)
  assign(msd, msdAux.dt)
  }


#Analisis parámetros ALL_devices

Parameters.dt <- fread(file.path(path_postLoc, paste0('Parameters_', parameters.method[1,2], '.csv')))

ggplot(Parameters.dt, aes(x=theta1, y=theta2)) + geom_point()


#Analisis MSDMASTER por experimento

comparacion <- msdmaster.dt %>% filter(experiment %in% c('4_cluster','5_cluster', 'All_devices'))
comparacion <- msdmaster.dt %>% filter(experiment %in% c('D_10%_4_cluster','D_10%_5_cluster', 'All_devices'))
comparacion <- msdmaster.dt %>% filter(experiment %in% c('D_25%_4_cluster','D_25%_5_cluster', 'All_devices'))
ggplot(comparacion, aes(x = experiment , y = cp_tp / tile_size)) +
  geom_boxplot(aes(fill = experiment)) +
  theme_bw() +
  labs(x = '', y = TeX('$b_{dt}$ (no. tiles)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

ggplot(comparacion, aes(x = time, y = cp_tp / tile_size, group = time)) +
  geom_boxplot(aes(fill = experiment)) +
  facet_grid(experiment ~ . ) +
  theme_bw() +
  labs(x = '\nTime (s)', y = TeX('$b_{dt}$ (no. tiles)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

ggplot(comparacion, aes(x = experiment, y = rd / tile_size)) +
  geom_boxplot(aes(fill = experiment)) +
  theme_bw() +
  labs(x = '', y = 'rmsd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


ggplot(comparacion, aes(x = time, y = rd / tile_size, group = time)) +
  geom_boxplot(aes(fill = experiment)) +
  facet_grid(experiment ~ . ) +
  theme_bw() +
  labs(x = '\nTime (s)', y = 'rmsd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


#Comparaciones por cluster de un experimento.

experimentcluster=`msd_D_10%_5_cluster`

ggplot(experimentcluster, aes(x = reorder(as.character(experimentcluster$cluster),experimentcluster$cluster) , y = cp_tp / tile_size)) +
  geom_boxplot(aes(fill = cluster)) +
  theme_bw() +
  labs(x = '', y = TeX('$b_{dt}$ (no. tiles)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

ggplot(experimentcluster, aes(x = time, y = cp_tp / tile_size, group = time)) +
  geom_boxplot(aes(fill = '')) +
  facet_grid(reorder(as.character(experimentcluster$cluster),experimentcluster$cluster) ~ . ) +
  theme_bw() +
  labs(x = '\nTime (s)', y = TeX('$b_{dt}$ (no. tiles)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

ggplot(experimentcluster, aes(x = reorder(as.character(experimentcluster$cluster),experimentcluster$cluster), y = rd / tile_size)) +
  geom_boxplot(aes(fill = cluster)) +
  theme_bw() +
  labs(x = '', y = 'rmsd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


ggplot(comparacion, aes(x = time, y = rd / tile_size, group = time)) +
  geom_boxplot(aes(fill = cluster)) +
  facet_grid(experiment ~ . ) +
  theme_bw() +
  labs(x = '\nTime (s)', y = 'rmsd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


#Histograma devices por cluster por experimento


devices_cluster = unique(experimentcluster[,cluster, by = 'device' ])

Parameterscluster = merge(Parameters.dt, devices_cluster, by = 'device')

ggplot(Parameterscluster, aes(x=theta1, y=theta2)) + 
  geom_point(aes(color = reorder(as.character(Parameterscluster$cluster),Parameterscluster$cluster))) +
  labs(colour = "Cluster") #+
     # geom_text(label=Parameterscluster$device)


hist(devices_cluster$cluster, main = 'Frecuencia de dispositivos por cluster',xlab = 'Cluster', ylab = 'Frecuencia', breaks=seq(min(devices_cluster$cluster)-0.5, max(devices_cluster$cluster)+0.5, by=1))

#Media de cambios de antena por cluster

devices_antenna_changes = unique(experimentcluster[,antenna_changes, by = c('device','cluster') ])

g = ggplot(devices_antenna_changes, aes(reorder(as.character(cluster),cluster), fill=reorder(as.character(antenna_changes-1),antenna_changes-1)) ) +
  labs(title = "Cambios de antena por cluster")+xlab("Cluster") + ylab("Cambios de antena")
  theme(plot.title = element_text(size = rel(2), colour = "blue"))

g+geom_bar(position="dodge") + 
  theme(axis.title.x = element_text(face="bold", size=10)) + 
  guides(fill=guide_legend(title="Cambios de antena"))

mean_antenna_changes =experimentcluster%>%
  group_by(cluster) %>%
  summarise(mean_antenna_changes = mean(antenna_changes))

barplot(as.matrix(mean_antenna_changes)[,2], main = 'Media de cambios de antena por cluster',xlab = "Cluster", ylab = "Media")

#Representacion antenas

ggplot(experimentcluster, aes(experimentcluster$event_cellID, experimentcluster$cp_tp)) +
  geom_point(aes(color = as.character(`msd_5_cluster`$cluster))) +
  labs(x = 'Eje x', y = 'Eje Y', title = 'Regresión', colour = 'Cluster') +
  theme_bw()



#Comparacion con regresion lineal entre dos experimentos deseados, mejor dejar que R autocomplete
#el experimento, ya que asi pone solo las comillas necesarias debido al %

experiment_X=msd_All_devices
experiment_Y=experimentcluster


comparingexperiments=merge(experiment_X,experiment_Y,by=c('device', 'time'))


lm.fit=lm(cp_tp.y ~ cp_tp.x, comparingexperiments)
summary(lm.fit)

ggplot(comparingexperiments, aes(comparingexperiments$cp_tp.x/tile_size, comparingexperiments$cp_tp.y/tile_size)) +
  geom_point(aes(color = as.character(comparingexperiments$cluster.y))) +
  labs(x = 'rmsd All devices', y = 'cp_tp 5 cluster', title = 'Regresión', colour = 'cluster') +
  stat_smooth(method = lm) +
  theme_bw() #+
   # geom_text(label=comparingexperiments$device)

confint(lm.fit)


lm.fit=lm(rd.y ~ rd.x, comparingexperiments)
summary(lm.fit)

ggplot(comparingexperiments, aes(comparingexperiments$rd.x/tile_size, comparingexperiments$rd.y/tile_size)) +
  geom_point(aes(color = as.character(comparingexperiments$cluster.y))) +
  labs(x = 'rmsd All devices', y = 'rmsd 5 cluster', title = 'Regresión', colour = 'cluster') +
  stat_smooth(method = lm) +
  theme_bw() #+
 # geom_text(label=comparingexperiments$device)

confint(lm.fit)


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####              PLOT DISTRIBUTIONS of b, rmsd AND msd                     ####
#### ** Bias                                                                ####
ggplot(msd_All_devices, aes(x = model, y = msd_All_devices$cp_tp / tile_size)) +
  geom_boxplot(aes(fill = model)) +
  theme_bw() +
  labs(x = '', y = TeX('$b_{dt}$ (no. tiles)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none') #+
  # geom_text(label=msd_All_devices$device)

ggplot(msd_All_devices, aes(x = time, y = msd_All_devices$cp_tp / tile_size, group = time)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nTime (s)', y = TeX('$b_{dt}$ (no. tiles)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

ggplot(msd.dt, aes(x = device, y = cp_tp / tile_size, group = device)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nDevice', y = TeX('$b_{dt}$ (no. tiles)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 6, angle = 90), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


#### ** rmsd 
ggplot(comparacion, aes(x = experiment, y = rd / tile_size)) +
  geom_boxplot(aes(fill = experiment)) +
  theme_bw() +
  labs(x = '', y = 'rmsd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


ggplot(comparacion, aes(x = time, y = rd / tile_size, group = time)) +
  geom_boxplot(aes(fill = experiment)) +
  facet_grid(experiment ~ . ) +
  theme_bw() +
  labs(x = '\nTime (s)', y = 'rmsd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

ggplot(msd.dt, aes(x = device, y = rd / tile_size, group = device)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nDevice', y = 'rmsd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 6, angle = 90), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


#### ** msd 
ggplot(msd.dt, aes(x = model, y = msd / tile_size**2)) +
  geom_boxplot(aes(fill = model)) +
  theme_bw() +
  labs(x = '', y = 'msd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


ggplot(msd.dt, aes(x = time, y = msd / tile_size**2, group = time)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nTime (s)', y = 'msd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


ggplot(msd.dt, aes(x = device, y = msd / tile_size**2, group = device)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nDevice', y = TeX('msd (no. tiles$^{2}$)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 6, angle = 90), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####            PLOT DEVICE TIME SERIES OF of b, rmsd AND msd               ####
devID <- sort(unique(msdmaster.dt$device))[30]
ggplot(msdmaster.dt[device == devID], aes(x = time, y = cp_tp / tile_size, group = time)) +
  geom_point(aes(color = experiment)) +
  facet_grid(model ~ . ) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  labs(x = '\nTime (s)', y = TeX('$b_{dt}$ (no. tiles)\n'), title = paste0('Device ', devID)) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12))
        # legend.position = 'none')

devID <- sort(unique(msd.dt$device))[5]
ggplot(msd.dt[device == devID], aes(x = time, y = rd / tile_size, group = time)) +
  geom_point(aes(color = model)) +
  facet_grid(model ~ . ) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  labs(x = '\nTime (s)', y = TeX('$b_{dt}$ (no. tiles)\n'), title = paste0('Device ', devID)) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

devID <- sort(unique(msd.dt$device))[5]
ggplot(msd.dt[device == devID], aes(x = time, y = msd / tile_size**2, group = time)) +
  geom_point(aes(color = model)) +
  facet_grid(model ~ . ) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  labs(x = '\nTime (s)', y = TeX('$b_{dt}$ (no. tiles)\n'), title = paste0('Device ', devID)) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')
