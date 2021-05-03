## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
dupProb <- read.csv(file = system.file('extdata/duplicity.csv', package = 'aggregation'))
head(dupProb)

## -----------------------------------------------------------------------------
regions <- read.csv(file = system.file('extdata/regions.csv', package = 'aggregation'))
head(regions)

## ----eval = FALSE-------------------------------------------------------------
#  # set the folder where the necessary input files are stored
#  path      <- 'extdata'
#  
#  prefix = 'postLocDevice'
#  
#  # set the duplicity probabilities file name, i.e. the file with duplicity probability for each device
#  dpFile <- system.file(path, 'duplicity.csv', package = 'aggregation')
#  
#  # set the regions file name, i.e. the file defining the regions for wich we need the estimation of the number
#  # of individuals detected by network.
#  rgFile <- system.file(path, 'regions.csv', package = 'aggregation')
#  
#  # set the path to the posterior location probabilities file
#  pathLoc <- system.file(path, package = 'aggregation')
#  
#  # set the number of random values to be generated
#  n <- 1e3
#  # call rNnetEvent
#  nNet <- rNnetEvent(n, dpFile, rgFile, pathLoc, prefix)
#  
#  head(nNet)

## ----eval = FALSE-------------------------------------------------------------
#  # print the mean number of detected individuals for each region, for each time instant
#  regions <- as.numeric(unique(nNet$region))
#  times <- unique(nNet$time)
#  
#  for(r in regions) {
#      print(paste0("region: ", r))
#      for(t in times) {
#          print(paste0("time instant: ", t, " number of individuals: " , mean(nNet[region == r][time ==t]$N)))
#      }
#  }

## ----eval = FALSE-------------------------------------------------------------
#  # For the origin-destination matrix we proceed as follows
#  # set the folder where the necessary input files are stored
#  path      <- 'extdata'
#  
#  prefixJ <- 'postLocJointProbDevice'
#  
#  # set the duplicity probabilities file name, i.e. the file with duplicity probability for each device
#  dpFile<-system.file(path, 'duplicity.csv', package = 'aggregation')
#  
#  # set the regions file name, i.e. the file defining the regions for wich we need the estimation of the number
#  # of individuals detected by network.
#  rgFile<-system.file(path, 'regions.csv', package = 'aggregation')
#  
#  # generate n random values
#  n <- 1e3
#  
#  nnetOD <- rNnetEventOD(n, dpFile, rgFile, system.file(path, package = 'aggregation'), prefixJ))
#  
#  head(nnetOD)

## ----eval = FALSE-------------------------------------------------------------
#  t_from <-0
#  t_to <- 10
#  
#  
#  regions_from <- sort(as.numeric(unique(nnetOD$region_from)))
#  regions_to <- sort(as.numeric(unique(nnetOD$region_to)))
#  
#  
#  ODmat <- matrix(nrow = length(regions_from), ncol = length(regions_to))
#  for(r1 in regions_from) {
#      for(r2 in regions_to) {
#            ODmat[r1,r2] <- round(mean(nnetOD[time_from==t1][time_to==t2][region_from==r1][region_to==r2]$Nnet))
#      }
#  }
#  ODmat

