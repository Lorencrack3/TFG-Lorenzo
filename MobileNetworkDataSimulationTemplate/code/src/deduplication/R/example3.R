#' Example of using deduplication package - showing the case of incompatibility between the type of connection in the
#' network and the emission model.
#'
#' This is an example of how to force the use of an emission model that differs from the type of connection done by the
#' network. Thanks to the information from the simulator and the flexibility of the functions of this package, it is
#' possible to make the assumption of a kind of emission model different from the real type of connection. Obviously,
#' this incoherence has consequences. It could cause to type of problems. Firstly there are events to antennas without 
#' signal under the assumption made to build the emission probability matrix. Secondly, there are transitions between 
#' pairs of antennas which are not possible in consecutive times because there is some
#' tile without signal between those antennas under the assumption made to build the emission probability matrix.  These
#' two problems can be detected by two functions in this package: \code{checkConnections_step1()} and
#' \code{checkConnections_step2()}. In the following example, one of these kinds of cases is shown.
#'
#' @examples
#'
#'
#' ###    READ DATA    ####
#'
#' # Set the folder where the necessary input files are stored
#'
#' path_root <- 'extdata'
#'
#'  # 0. Read simulation params
#'
#' simParams <- readSimulationParams(system.file(path_root, 'simulation.xml', package = 'deduplication'))
#'
#'  # 1. Read grid parameters
#'
#' gridParams <- readGridParams(system.file(path_root, 'grid.csv', package = 'deduplication'))
#'
#'  # 2. Read network events
#'
#' events <- readEvents(system.file(path_root, 'AntennaInfo_MNO_MNO1.csv', package = 'deduplication'))
#'
#' # 3. Set the signal file name, i.e. the file where the signal strength/quality for each tile in the grid is stored
#'
#' signalFileName <- system.file(path_root, 'SignalMeasure_MNO1.csv')
#'
#'  # 4. Set the antennas file name
#'
#' antennasFileName <- system.file(path_root, 'antennas.xml')
#'
#' ####    PREPARE DATA    ####
#'
#'  # 5. Initial state distribution (prior)
#'
#' nTiles <- gridParams$ncol * gridParams$nrow initialDistr_RSS_uniform.vec <- initialDistr_SDM_uniform.vec <- rep(1 /
#' nTiles, nTiles)
#'
#'  # 6. Get a list of detected devices
#'
#'  devices <- getDeviceIDs(events)
#'
#'  # 7. Get connections for each device
#'
#' connections <- getConnections(events)
#'
#' # 8. Emission probabilities are computed from the signal strength/quality file.  In this case handoverType is strength then the emission model should be RSS.
#'
#' emissionProbs_RSS <- getEmissionProbs(nrows = gridParams$nrow, ncols = gridParams$ncol, signalFileName = signalFileName, sigMin = simParams$conn_threshold, handoverType = simParams$connection_type[[1]], emissionModel = "RSS", antennaFileName = antennasFileName)
#'
#'  # We also try to use the emission model with the SDM assumption.
#'
#' emissionProbs_SDM <- getEmissionProbs(nrows = gridParams$nrow, ncols = gridParams$ncol, signalFileName = signalFileName, sigMin = simParams$conn_threshold, handoverType = simParams$connection_type[[1]], emissionModel = "SDM", antennaFileName = antennasFileName)
#'
#'  # 9. Build joint emission probabilities for both options: RSS and SDM.
#'
#' jointEmissionProbs_RSS <- getEmissionProbsJointModel(emissionProbs_RSS)
#'
#' jointEmissionProbs_SDM <- getEmissionProbsJointModel(emissionProbs_SDM)
#'
#'  # 10. Build the generic model for both cases and the a priori uniform (by default)
#'
#' model_RSS_uniform <- getGenericModel( nrows = gridParams$nrow, ncols = gridParams$ncol, emissionProbs_RSS, initSteady = FALSE, aprioriProb = initialDistr_RSS_uniform.vec)
#'
#' model_SDM_uniform <- getGenericModel( nrows = gridParams$nrow, ncols = gridParams$ncol, emissionProbs_SDM, initSteady = FALSE, aprioriProb = initialDistr_SDM_uniform.vec)
#'
#'  # 11. Fit models.
#'
#' ll_RSS_uniform <- fitModels(length(devices), model_RSS_uniform, connections)
#'
#' ll_SDM_uniform <- fitModels(length(devices), model_SDM_uniform, connections)
#'
#' # The log likelihood is infinity for those connections that are impossible under the model SDM (ll_SDM_uniform).
#' 
#'
#'  # 12. Make the checks and corresponding imputations to force the fit of SDM model.
#'
#' ### check1: connections incompatibles with emissionProbs are imputed as NA
#'
#' check1_SDM <- checkConnections_step1(connections, emissionProbs_SDM)
#'
#' check1_SDM$infoCheck_step1 connections_ImpSDM0 <-check1_SDM$connectionsImp
#'
#' ### check2: jumps incompatibles with emissionProbs are avoid by making time padding.
#'
#' check2_SDM <- checkConnections_step2(emissionProbs = emissionProbs_SDM, connections = connections_ImpSDM0, gridParams = gridParams)
#'
#' connections_ImpSDM <- check2_SDM$connections_pad
#'
#'
#'  # 13. Fit models.
#'
#' ll_SDM_uniform <- fitModels(length(devices), model_SDM_uniform, connections_ImpSDM)
#'
#'  # 14. Build the joint model.
#'
#' modelJ_RSS_uniform <- getJointModel( nrows = gridParams$nrow, ncols = gridParams$ncol, jointEmissionProbs = jointEmissionProbs_RSS, initSteady = FALSE, aprioriJointProb = initialDistr_RSS_uniform.vec)
#'
#' modelJ_SDM_uniform <- getJointModel( nrows = gridParams$nrow, ncols = gridParams$ncol, jointEmissionProbs = jointEmissionProbs_SDM, initSteady = FALSE, aprioriJointProb = initialDistr_SDM_uniform.vec)
#'
#'  # 15. Apriori probability of 2-to-1
#'
#' Pii <- aprioriOneDeviceProb(simParams$prob_sec_mobile_phone, length(devices))
#'
#'  # 16. Build a matrix of pairs of devices to compute duplicity probability
#'
#' pairs4dup_RSS <- computePairs(connections, length(devices), oneToOne = TRUE)
#'
#' pairs4dup_SDM <- computePairs(connections_ImpSDM, length(devices), oneToOne = TRUE)
#'
#'
#' ####    COMPUTE DUPLICITY    ####
#'
#'
#' # 17. Compute duplicity probabilities using "1to1" method
#'
#' duplicity1to1_RSS_uniform.dt <- computeDuplicityBayesian("1to1", devices, pairs4dup_RSS, modelJ_RSS_uniform,
#' ll_RSS_uniform, P1 = NULL, Pii=Pii)
#'
#' duplicity1to1_SDM_uniform.dt <- computeDuplicityBayesian("1to1", devices,
#' pairs4dup_SDM, modelJ_SDM_uniform, ll_SDM_uniform, P1 = NULL, Pii=Pii)
example3 <- function() {}
