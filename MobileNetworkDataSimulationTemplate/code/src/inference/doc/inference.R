## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  library(inference, warn.conflicts = FALSE)
#  
#  path      <- 'extdata'
#  
#  prefix <- 'postLocDevice'
#  postLocPath <- system.file(path, package = 'inference')
#  
#  
#  dpFileName <- system.file(path, 'duplicity.csv', package = 'inference')
#  rgFileName <- system.file(path, 'regions.csv', package = 'inference')
#  
#  omega_r <- computeDeduplicationFactors(dupFileName = dpFileName,
#                                         regsFileName = rgFileName,
#                                         postLocPrefix = prefix,
#                                         postLocPath = postLocPath)
#  head(omega_r)

## ----eval = FALSE-------------------------------------------------------------
#  pRFileName <- system.file(path, 'pop_reg.csv', package = 'inference')
#  pRateFileName <- system.file(path, 'pnt_rate.csv', package = 'inference')
#  grFileName <- system.file(path, 'grid.csv', package = 'inference')
#  params <- computeDistrParams(omega = omega_r,
#                               popRegFileName = pRFileName,
#                               pntRateFileName = pRateFileName,
#                               regsFileName = rgFileName,
#                               gridFileName = grFileName)
#  head(params)

## ----eval = FALSE-------------------------------------------------------------
#  nFileName <- system.file(path, 'nnet.csv', package = 'inference')
#  nnet <- readNnetInitial(nFileName)
#  
#  # Beta Negative Binomial distribution
#  n_bnb <- computeInitialPopulation(nnet = nnet,
#                                    params = params,
#                                    popDistr = 'BetaNegBin',
#                                    rndVal = TRUE)
#  
#  head(n_bnb$stats)
#  head(n_bnb$rnd_values)

## ----eval = FALSE-------------------------------------------------------------
#  # Negative Binomial distribution
#  n_nb <- computeInitialPopulation(nnet = nnet,
#                                   params = params,
#                                   popDistr = 'NegBin',
#                                   rndVal = TRUE)
#  
#  head(n_nb$stats)
#  head(n_nb$rnd_values)

## ----eval = FALSE-------------------------------------------------------------
#  # State process Negative Binomial distribution
#  n_stnb <- computeInitialPopulation(nnet = nnet,
#                                    params = params,
#                                    popDistr= 'STNegBin',
#                                    rndVal = TRUE)
#  
#  head(n_stnb$stats)
#  head(n_stnb$rnd_values)
#  

## ----eval = FALSE-------------------------------------------------------------
#  nnetODFile <- system.file(path, 'nnetOD.zip', package = 'inference')

## ----eval=FALSE---------------------------------------------------------------
#  # Beta Negative Binomial distribution
#  nt_bnb <- computePopulationT(nt0 = n_bnb$rnd_values,
#                               nnetODFileName = nnetODFile,
#                               rndVal = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  
#  times <- names(nt_bnb)
#  t <- sample(1:length(times), size = 1)
#  t
#  head(nt_bnb[[t]]$stats)
#  head(nt_bnb[[t]]$rnd_values)

## ----eval=FALSE---------------------------------------------------------------
#  # Negative Binomial distribution
#  nt_nb <- computePopulationT(nt0 = n_nb$rnd_values,
#                              nnetODFileName = nnetODFile,
#                              rndVal = TRUE)
#  
#  # to display results, select a random time instant
#  times <- names(nt_nb)
#  t <- sample(1:length(times), size = 1)
#  t
#  head(nt_nb[[t]]$stats)
#  head(nt_nb[[t]]$rnd_values)

## ----eval=FALSE---------------------------------------------------------------
#  # State process Negative Binomial distribution
#  nt_stnb <- computePopulationT(nt0 = n_stnb$rnd_values,
#                                nnetODFileName = nnetODFile,
#                                rndVal = TRUE)
#  
#  # to display results, select a random time instant
#  times <- names(nt_stnb)
#  t <- sample(1:length(times), size = 1)
#  t
#  head(nt_stnb[[t]]$stats)
#  head(nt_stnb[[t]]$rnd_values)

## ----eval=FALSE, message=FALSE------------------------------------------------
#  # Beta Negative Binomial distribution
#  OD_bnb <- computePopulationOD(nt0 = n_bnb$rnd_values,
#                                nnetODFileName = nnetODFile,
#                                rndVal = TRUE)
#  
#  # to display results, select a random time instant
#  time_pairs <- names(OD_bnb)
#  i <- sample(1:length(time_pairs), size = 1)
#  time_pairs[i]
#  head(OD_bnb[[i]]$stats)
#  head(OD_bnb[[i]]$rnd_values)

## ----eval=FALSE---------------------------------------------------------------
#  # Negative Binomial distribution
#  OD_nb <- computePopulationOD(nt0 = n_nb$rnd_values,
#                               nnetODFileName = nnetODFile,
#                               rndVal = TRUE)
#  
#  # to display results, select a random time instant
#  time_pairs <- names(OD_nb)
#  i <- sample(1:length(time_pairs), size = 1)
#  time_pairs[i]
#  head(OD_nb[[i]]$stats)
#  head(OD_nb[[i]]$rnd_values)

## ----eval=FALSE---------------------------------------------------------------
#  # State process Negative Binomial distribution
#  OD_stnb <- computePopulationOD(nt0 = n_stnb$rnd_values,
#                                 nnetODFileName = nnetODFile,
#                                 rndVal = TRUE)
#  
#  # to display results, select a random time instant
#  time_pairs <- names(OD_stnb)
#  i <- sample(1:length(time_pairs), size = 1)
#  time_pairs[i]
#  head(OD_stnb[[i]]$stats)
#  head(OD_stnb[[i]]$rnd_values)
#  

## ----eval=FALSE, message=FALSE------------------------------------------------
#  library(plumber)
#  pathPL <-'plumber'
#  initPop_api <- plumber::plumb(system.file(pathPL, 'plumb.R', package = 'inference'))
#  initPop_api$run(host = '127.0.0.1', port = 8000)

## ----eval=FALSE, message=FALSE------------------------------------------------
#  library(httr)
#  library(jsonlite)
#  library(data.table)
#  
#  # set the folder where the necessary input files are stored and the prefix of the input file names.
#  path      <- 'extData'
#  prefix <- 'postLocDevice'
#  dpFileName <- system.file(path, 'duplicity.csv', package = 'inference')
#  rgFileName <- system.file(path, 'regions.csv', package = 'inference')
#  

## ----eval=FALSE, message=FALSE------------------------------------------------
#  # prepare the body of the http request
#  body <- list(
#   .dupFileName = dpFileName,
#   .regsFileName = rgFileName,
#   .postLocPrefix = prefix,
#   .postLocPath = postLocPath
#  )
#  
#  # set API path
#  pathDedup <- 'computeDeduplicationFactors'
#  
#  # send POST Request to API
#  url <- "http://127.0.0.1:8000"
#  raw.result <- POST(url = url, path = pathDedup, body = body, encode = 'json')

## ----eval = FALSE-------------------------------------------------------------
#  # check status code
#  raw.result$status_code
#  
#  # transform back the results from json format
#  omega_r <- as.data.table(fromJSON(rawToChar(raw.result$content)))

## ----eval = FALSE-------------------------------------------------------------
#  # Compute the parameters of the distribution
#  # First reads the number of individuals detected by network
#  nFileName <- system.file(path, 'nnet.csv', package = 'inference')
#  nnet <- readNnetInitial(nFileName)
#  
#  pRFileName <- system.file(path, 'pop_reg.csv', package = 'inference')
#  pRateFileName <- system.file(path, 'pnt_rate.csv', package = 'inference')
#  grFileName <- system.file(path, 'grid.csv', package = 'inference')
#  
#  # prepare the body of the http request
#  body <- list(
#    .omega = omega_r,
#   .popRegFileName = pRFileName,
#   .pntRateFileName  = pRateFileName,
#   .regsFileName = rgFileName,
#   .gridFileName = grFileName,
#   .rel_bias = 0,
#   .cv = 1e-5
#  )
#  
#  # set API path
#  pathDistr <- 'computeDistrParams'
#  
#  # send POST Request to API
#  raw.result <- POST(url = url, path = pathDistr, body = body, encode = 'json')
#  
#  # check status code
#  raw.result$status_code
#  
#  # transform back the results from json format
#  params <- as.data.table(fromJSON(rawToChar(raw.result$content)))

## ----eval = FALSE-------------------------------------------------------------
#  # Compute the population count distribution at t0 using the Beta Negative Binomial distribution
#  
#  # prepare the body of the http request
#  body <- list(
#   .nnet = nnet,
#   .params = params,
#   .popDistr = 'BetaNegBin',
#   .rndVal = TRUE,
#   .ciprob = 0.95,
#   .method = 'ETI'
#  )
#  
#  # set API path
#  pathInit <- 'computeInitialPopulation'
#  
#  # send POST Request to API
#  raw.result <- POST(url = url, path = pathInit, body = body, encode = 'json')
#  
#  # check status code
#  raw.result$status_code
#  
#  # transform back the results from json format
#  n_bnb <- fromJSON(rawToChar(raw.result$content))
#  
#  # display results
#  n_nb$stats
#  head(n_nb$rnd_values)

## ----eval = FALSE-------------------------------------------------------------
#  # Compute the population count distribution at time instants t > t0 using the Beta Negative Binomial distribution
#  # first set the name of the file with the population moving from one region to another (output of the aggregation # package)
#  nnetODFile <- system.file(path, 'nnetOD.zip', package = 'inference')
#  
#  # prepare the body of the http request
#  body <- list(
#   .nt0 = as.data.table(n_bnb$rnd_values),
#   .nnetODFileName = nnetODFile,
#   .zip = TRUE,
#   .rndVal = TRUE,
#   .ciprob = 0.95,
#   .method = 'ETI'
#  )
#  
#  # set API path
#  pathT <- 'computePopulationT'
#  
#  # send POST Request to API
#  raw.result <- POST(url = url, path = pathT, body = body, encode = 'json')
#  
#  # check status code
#  raw.result$status_code
#  
#  # transform back the results from json format
#  nt_bnb <- fromJSON(rawToChar(raw.result$content))
#  
#  # display results
#  # first, select a random time instant
#  times <- names(nt_bnb)
#  t <- sample(1:length(times), size = 1)
#  t
#  nt_bnb[[t]]$stats
#  head(nt_bnb[[t]]$rnd_values)

## ----eval = FALSE-------------------------------------------------------------
#  # Compute the Origin-Destination matrices for all pairs of time instants time_from-time_to using the Beta Negative # Binomial distribution
#  
#  # prepare the body of the http request
#  body <- list(
#   .nt0 = as.data.table(n_bnb$rnd_values),
#   .nnetODFileName = nnetODFile,
#   .zip = TRUE,
#   .rndVal = TRUE,
#   .ciprob = 0.95,
#   .method = 'ETI'
#  )
#  
#  # set API path
#  pathOD <- 'computePopulationOD'
#  
#  # send POST Request to API
#  raw.result <- POST(url = url, path = pathOD, body = body, encode = 'json')
#  
#  # check status code
#  raw.result$status_code
#  
#  # transform back the results from json format
#  OD_bnb <- fromJSON(rawToChar(raw.result$content))
#  
#  # display results
#  time_pairs <- names(OD_bnb)
#  # first, select a random time instants pair
#  i <- sample(1:length(time_pairs), size = 1)
#  time_pairs[i]
#  OD_bnb[[i]]$stats
#  head(OD_bnb[[i]]$rnd_values)

