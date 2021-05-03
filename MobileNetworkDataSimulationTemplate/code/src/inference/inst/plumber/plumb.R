
#' @post /computeInitialPopulation
#' @get /computeInitialPopulation
function(req) {
    require(jsonlite)
    require(inference)
    require(data.table)
    # post body
    body <- jsonlite::fromJSON(req$postBody)
    
    .nnet <- as.data.table(body$.nnet)
    .params <- as.data.table(body$.params)
    .popDistr <- body$.popDistr
    .rndVal <- body$.rndVal
    .ciprob <- body$.ciprob
    .method <- body$.method
    
    inference::computeInitialPopulation(nnet = .nnet, params = .params, popDistr = .popDistr, rndVal = .rndVal, ciprob = .ciprob, method = .method)
}

#' @post /computePopulationT
#' @get /computePopulationT
function(req) {
    require(jsonlite)
    require(inference)
    require(data.table)
    # post body
    body <- jsonlite::fromJSON(req$postBody)
    
    .nt0 <- as.data.table(body$.nt0)
    .nnetODFileName <- body$.nnetODFileName
    .zip <- body$.zip
    .rndVal <- body$.rndVal
    .ciprob <- body$.ciprob
    .method <- body$.method
    
    inference::computePopulationT(nt0 = .nt0, nnetODFileName = .nnetODFileName, zip = .zip, rndVal = .rndVal, ciprob = .ciprob, method = .method)
}

#' @post /computePopulationOD
#' @get /computePopulationOD
function(req) {
    require(jsonlite)
    require(inference)
    require(data.table)
    # post body
    body <- jsonlite::fromJSON(req$postBody)
    
    .nt0 <- as.data.table(body$.nt0)
    .nnetODFileName <- body$.nnetODFileName
    .zip <-body$.zip
    .rndVal <- body$.rndVal
    .ciprob <- body$.ciprob
    .method <- body$.method
    
    inference::computePopulationOD(nt0 = .nt0, nnetODFileName = .nnetODFileName, zip = .zip, rndVal = .rndVal, ciprob = .ciprob, method = .method)
}

#' @post /computeDeduplicationFactors
#' @get /computeDeduplicationFactors
function(req) {
    require(jsonlite)
    require(inference)
    require(data.table)
    # post body
    body <- jsonlite::fromJSON(req$postBody)
    
    .dupFileName <- body$.dupFileName
    .regsFileName <- body$.regsFileName
    .postLocPrefix <-body$.postLocPrefix
    .postLocPath <- body$.postLocPath

    inference::computeDeduplicationFactors(dupFileName = .dupFileName, regsFileName = .regsFileName, postLocPrefix = .postLocPrefix, postLocPath = .postLocPath) 
}

#' @post /computeDistrParams
#' @get /computeDistrParams
function(req) {
    require(jsonlite)
    require(inference)
    require(data.table)
    # post body
    body <- jsonlite::fromJSON(req$postBody)
    
    .omega <- as.data.table(body$.omega)
    .popRegFileName <- body$.popRegFileName
    .pntRateFileName <-body$.pntRateFileName
    .regsFileName <- body$.regsFileName
    .gridFileName <- body$.gridFileName
    .rel_bias <-body$.rel_bias
    .cv <- body$.cv
    
    inference::computeDistrParams(omega = .omega, popRegFileName = .popRegFileName, pntRateFileName = .pntRateFileName, regsFileName = .regsFileName, gridFileName = .gridFileName, rel_bias = .rel_bias, cv = .cv)

}

#' @post /getPath
#' @get /getPath
function(req) {
    require(jsonlite)
    # post body
    path <- jsonlite::fromJSON(req$postBody)
    system.file(path, package = 'inference')
    
}

#' @post /getNnetInitial
#' @get /getNnetInitial
function(req) {
    require(jsonlite)
    # post body
    body <- jsonlite::fromJSON(req$postBody)
    readNnetInitial(body$.nnetFileName)
}
