# ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE--------------------------------------------------------------
library(deduplication)
path_root <- 'extdata'

## -----------------------------------------------------------------------------
gridParams <-readGridParams(system.file(path_root, 'grid.csv', package = 'deduplication'))```


## -----------------------------------------------------------------------------
simParams <-readSimulationParams(system.file(path_root, 'simulation.xml', package = 'deduplication'))

## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
devices <- getDeviceIDs(events)
connections <- getConnections(events)

## -----------------------------------------------------------------------------
emissionProbs <- getEmissionProbs(gridParams$nrow, gridParams$ncol, system.file(path_root, 'SignalMeasure_MNO1.csv', package = 'deduplication'), simParams$conn_threshold)
jointEmissionProbs <- getEmissionProbsJointModel(emissionProbs)

## -----------------------------------------------------------------------------
model <- getGenericModel(gridParams$nrow, gridParams$ncol, emissionProbs)
modelJ <- getJointModel(gridParams$nrow, gridParams$ncol, jointEmissionProbs)

## -----------------------------------------------------------------------------
ll <- fitModels(length(devices), model, connections)

## -----------------------------------------------------------------------------
coverarea <- readCells(system.file(path_root, 'AntennaCells_MNO1.csv', package = 'deduplication'))
antennaNeigh <- antennaNeighbours(coverarea)

## -----------------------------------------------------------------------------
P1 <- aprioriDuplicityProb(simParams$prob_sec_mobile_phone, length(devices))

## -----------------------------------------------------------------------------
pairs4dup<-computePairs(connections, length(devices), oneToOne = FALSE, P1 = P1, limit = 0.05, antennaNeighbors = antennaNeigh)

## -----------------------------------------------------------------------------
probDup <- computeDuplicityBayesian("pairs", devices, pairs4dup, modelJ, ll, P1)

## -----------------------------------------------------------------------------
gridParams <-readGridParams(system.file(path_root, 'grid.csv', package = 'deduplication'))
simParams <-readSimulationParams(system.file(path_root, 'simulation.xml', package = 'deduplication'))
events <- readEvents(system.file(path_root, 'AntennaInfo_MNO_MNO1.csv', package = 'deduplication'))
devices <- getDeviceIDs(events)
connections <- getConnections(events)
emissionProbs <- getEmissionProbs(gridParams$nrow, gridParams$ncol, system.file(path_root, 'SignalMeasure_MNO1.csv', package = 'deduplication'), simParams$conn_threshold)
jointEmissionProbs <- getEmissionProbsJointModel(emissionProbs)
model <- getGenericModel(gridParams$nrow, gridParams$ncol, emissionProbs)
modelJ <- getJointModel(gridParams$nrow, gridParams$ncol, jointEmissionProbs)
ll <- fitModels(length(devices), model, connections)

## -----------------------------------------------------------------------------
Pii <- aprioriOneDeviceProb(simParams$prob_sec_mobile_phone, length(devices))
pairs4dup<-computePairs(connections, length(devices), oneToOne = TRUE)

## -----------------------------------------------------------------------------
probDup2 <- computeDuplicityBayesian("1to1", devices, pairs4dup, modelJ, ll, P1 = NULL, Pii=Pii)

## -----------------------------------------------------------------------------
probDup3 <- computeDuplicityBayesian(method, devices, pairs4dup, modelJ, ll, P1 = NULL, Pii = NULL, init = TRUE, lambda = 0.67)

## -----------------------------------------------------------------------------
gridParams <-readGridParams(system.file(path_root, 'grid.csv', package = 'deduplication'))
events <- readEvents(system.file(path_root, 'AntennaInfo_MNO_MNO1.csv', package = 'deduplication'))
devices <- getDeviceIDs(events)
T<-nrow(unique(events[,1]))
coverarea <- readCells(system.file(path_root, 'AntennaCells_MNO1.csv', package = 'deduplication'))
antennaNeigh <- antennaNeighbours(coverarea)
P1a <- aprioriDuplicityProb(simParams$prob_sec_mobile_phone, length(devices))
pairs4dup<-computePairs(connections, length(devices), oneToOne = FALSE, antennaNeighbors = antennaNeigh)
probDup3 <-computeDuplicityTrajectory(path=path_root, devices, gridParams, pairs4dup, P1 = P1a, T, gamma = 0.5)


## -----------------------------------------------------------------------------
path_root <- 'extdata'

## -----------------------------------------------------------------------------
gridfile <- system.file(path_root, 'grid.csv', package = 'deduplication')

## -----------------------------------------------------------------------------
eventsfile <- system.file(path_root, 'AntennaInfo_MNO_MNO1.csv', package = 'deduplication')

## -----------------------------------------------------------------------------
signalfile<-system.file(path_root, 'SignalMeasure_MNO1.csv', package = 'deduplication')

## -----------------------------------------------------------------------------
antennacellsfile <- system.file(path_root, 'AntennaCells_MNO1.csv', package = 'deduplication')

## -----------------------------------------------------------------------------
simulationfile<-system.file(path_root, 'simulation.xml', package = 'deduplication')
```

## -----------------------------------------------------------------------------
out1<-computeDuplicity("pairs", gridFileName = gridfile, eventsFileName = eventsfile, signalFileName = signalfile, antennaCellsFileName = antennacellsfile, simulationFileName = simulationfile)

## -----------------------------------------------------------------------------
out2<-computeDuplicity("1to1", gridFileName = gridfile, eventsFileName = eventsfile, signalFileName = signalfile, simulatedData = TRUE, simulationFileName = simulationfile)

## -----------------------------------------------------------------------------
out2p<-computeDuplicity("1to1", gridFileName = gridfile, eventsFileName = eventsfile, signalFileName = signalfile, simulatedData = TRUE, simulationFileName = simulationfile, lambda = 0.67)

## -----------------------------------------------------------------------------
out3<-computeDuplicity("trajectory", gridFileName = gridfile, eventsFileName = eventsfile, signalFileName = signalfile, antennaCellsFileName = antennacellsfile, simulationFileName = simulationfile, path= path_root)

## -----------------------------------------------------------------------------
aprioriProbModel <- matrix (1 / (gridParams$nrow * gridParams$ncol), nrow = gridParams$nrow, ncol = gridParams$ncol)
model <- getGenericModel(gridParams$nrow, gridParams$ncol, emissionProbs, initSteady = FALSE, aprioriProb = aprioriProbModel)
modelJ <- getJointModel(gridParams$nrow, gridParams$ncol, jointEmissionProbs, initSteady = FALSE, aprioriJointProb = aprioriProbModel)
