#' Example of using the API inference
#'
#' # This is an example on how to use the REST API of this package to generate the distribution of the population count.
#'
#' @examples
#'
#' ########################################
#' # First, in a separate R console run the following instructions that start the http API on your local computer
#' # One can replace 127.0.0.1 with the API address of another Web server
#' library(plumber)
#' pathPL <-'plumber'
#' initPop_api <- plumber::plumb(system.file(pathPL, 'plumb.R', package = 'inference'))
#' initPop_api$run(host = '127.0.0.1', port = 8000)
#' ########################################
#'
#' library(httr)
#' library(jsonlite)
#' library(data.table)
#'
#' # set the folder where the necessary input files are stored and the prefix of the input file names.
#' path      <- 'extData'
#' prefix <- 'postLocDevice'
#' postLocPath <- system.file(path, package = 'inference')
#'
#' # compute the deduplication factors
#' dpFileName <- system.file(path, 'duplicity.csv', package = 'inference')
#' rgFileName <- system.file(path, 'regions.csv', package = 'inference')
#'
#' # prepare the body of the http request
#' body <- list(
#'   .dupFileName = dpFileName,
#'   .regsFileName = rgFileName,
#'   .postLocPrefix = prefix,
#'   .postLocPath = postLocPath
#' )
#'
#' # set API path
#' pathDedup <- 'computeDeduplicationFactors'
#'
#' # send POST Request to API
#' url <- "http://127.0.0.1:8000"
#' raw.result <- POST(url = url, path = pathDedup, body = body, encode = 'json')
#'
#' # check status code
#' raw.result$status_code
#'
#' # transform back the results from json format
#' omega_r <- as.data.table(fromJSON(rawToChar(raw.result$content)))
#'
#' # Compute the parameters of the distribution
#' # First reads the number of individuals detected by network
#' nFileName <- system.file(path, 'nnet.csv', package = 'inference')
#' nnet <- readNnetInitial(nFileName)
#'
#' pRFileName <- system.file(path, 'pop_reg.csv', package = 'inference')
#' pRateFileName <- system.file(path, 'pnt_rate.csv', package = 'inference')
#' grFileName <- system.file(path, 'grid.csv', package = 'inference')
#'
#' # prepare the body of the http request
#' body <- list(
#'    .omega = omega_r,
#'   .popRegFileName = pRFileName,
#'   .pntRateFileName  = pRateFileName,
#'   .regsFileName = rgFileName,
#'   .gridFileName = grFileName,
#'   .rel_bias = 0,
#'   .cv = 1e-5
#' )
#'
#' # set API path
#' pathDistr <- 'computeDistrParams'
#'
#'  # send POST Request to API
#' raw.result <- POST(url = url, path = pathDistr, body = body, encode = 'json')
#'
#' # check status code
#' raw.result$status_code
#'
#' # transform back the results from json format
#' params <- as.data.table(fromJSON(rawToChar(raw.result$content)))
#'
#' # Compute the population count distribution at t0 using the Beta Negative Binomial distribution
#'
#' # prepare the body of the http request
#' body <- list(
#'   .nnet = nnet,
#'   .params = params,
#'   .popDistr = 'BetaNegBin',
#'   .rndVal = TRUE,
#'   .ciprob = 0.95,
#'   .method = 'ETI'
#' )
#'
#' # set API path
#' pathInit <- 'computeInitialPopulation'
#'
#'  # send POST Request to API
#' raw.result <- POST(url = url, path = pathInit, body = body, encode = 'json')
#'
#' # check status code
#' raw.result$status_code
#'
#' # transform back the results from json format
#' n_bnb <- fromJSON(rawToChar(raw.result$content))
#'
#' # display results
#' n_bnb$stats
#' head(n_bnb$rnd_values)
#'
#' # Compute the population count distribution at time instants t > t0 using the Beta Negative Binomial distribution
#' # first set the name of the file with the population moving from one region to another (output of the aggregation package)
#' nnetODFile <- system.file(path, 'nnetOD.zip', package = 'inference')
#'
#' # prepare the body of the http request
#' body <- list(
#'   .nt0 = as.data.table(n_bnb$rnd_values),
#'   .nnetODFileName = nnetODFile,
#'   .zip = TRUE,
#'   .rndVal = TRUE,
#'   .ciprob = 0.95,
#'   .method = 'ETI'
#' )
#'
#' # set API path
#' pathT <- 'computePopulationT'
#'
#'  # send POST Request to API
#' raw.result <- POST(url = url, path = pathT, body = body, encode = 'json')
#'
#' # check status code
#' raw.result$status_code
#'
#' # transform back the results from json format
#' nt_bnb <- fromJSON(rawToChar(raw.result$content))
#'
#' # display results
#' # first, select a random time instant
#' times <- names(nt_bnb)
#' t <- sample(1:length(times), size = 1)
#' t
#' nt_bnb[[t]]$stats
#' head(nt_bnb[[t]]$rnd_values)
#'
#'
#' # Compute the Origin-Destination matrices for all pairs of time instants time_from-time_to using the Beta Negative Binomial distribution
#'
#' # prepare the body of the http request
#' body <- list(
#'   .nt0 = as.data.table(n_bnb$rnd_values),
#'   .nnetODFileName = nnetODFile,
#'   .zip = TRUE,
#'   .rndVal = TRUE,
#'   .ciprob = 0.95,
#'   .method = 'ETI'
#' )
#'
#' # set API path
#' pathOD <- 'computePopulationOD'
#'
#'  # send POST Request to API
#' raw.result <- POST(url = url, path = pathOD, body = body, encode = 'json')
#'
#' # check status code
#' raw.result$status_code
#'
#' # transform back the results from json format
#' OD_bnb <- fromJSON(rawToChar(raw.result$content))
#'
#' # display results
#' time_pairs <- names(OD_bnb)
#' # first, select a random time instants pair
#' i <- sample(1:length(time_pairs), size = 1)
#' time_pairs[i]
#' OD_bnb[[i]]$stats
#' head(OD_bnb[[i]]$rnd_values)

exampleAPI <- function() {}
