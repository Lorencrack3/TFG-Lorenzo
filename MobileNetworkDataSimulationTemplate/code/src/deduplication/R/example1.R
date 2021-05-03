#' Example of using deduplication package - the simple way
#' 
#' This is just an example on how to compute duplicity probabilities using simulated data. All the files used in this
#' example are supposed to be produced using the simulation software. The "simulation.xml" file is an exeception and it
#' is an input file for the simulation software. The files used in this example are provided with the deduplication
#' package.
#'
#' @examples
#'
#' # set the folder where the necessary input files are stored
#' path_root      <- 'extdata'
#' 
#' # set the grid file name, i.e. the file where the grid parameters are found
#' gridfile <-system.file(path_root, 'grid.csv', package = 'deduplication')
#' 
#' # set the events file name, i.e. the file with network events registered during a simulation
#' eventsfile<-system.file(path_root, 'AntennaInfo_MNO_MNO1.csv', package = 'deduplication')
#'
#' # set the signal file name, i.e. the file where the signal strength/quality for each tile in the grid is stored
#' signalfile<-system.file(path_root, 'SignalMeasure_MNO1.csv', package = 'deduplication')
#'
#' # set the antenna cells file name, i.e. the file where the simulation software stored the coverage area for each antenna
#' # This file is needed only if the duplicity probabilities are computed using "pairs" method
#' antennacellsfile<-system.file(path_root, 'AntennaCells_MNO1.csv', package = 'deduplication')
#'
#' # set the simulation file name, i.e. the file with the simulation parameters used to produce the data set
#' simulationfile<-system.file(path_root, 'simulation.xml', package = 'deduplication')
#' 
#' # compute the duplicity probabilities using the "pairs" method
#' out1<-computeDuplicity("pairs", gridFileName = gridfile, eventsFileName = eventsfile, signalFileName = signalfile, antennaCellsFileName = antennacellsfile, simulationFileName = simulationfile)
#'
#' # compute the duplicity probabilities using the "1to1" method
#' out2<-computeDuplicity("1to1", gridFileName = gridfile, eventsFileName = eventsfile, signalFileName = signalfile, simulatedData = TRUE, simulationFileName = simulationfile)
#' 
#' # compute the duplicity probabilities using the "1to1" method with lambda
#' out2p<-computeDuplicity("1to1", gridFileName = gridfile, eventsFileName = eventsfile, signalFileName = signalfile, simulatedData = TRUE, simulationFileName = simulationfile, lambda = 0.67)
#' 
#' # compute the duplicity probabilities using the "trajectory" method
#' prefix <- 'postLocDevice'
#' out3<-computeDuplicity("trajectory", gridFileName = gridfile, eventsFileName = eventsfile, signalFileName = signalfile, antennaCellsFileName = antennacellsfile, simulationFileName = simulationfile, path= system.file(path_root, package='deduplication'), prefix = prefix)
#'
example1 <- function() {}
