#' @title Reads the parameters of the simulation used to generate a data set.
#'
#' @description Reads the parameters of the simulation used to generate a data set from an .xml file used by the
#'   simulation software. The following parameters are needed by this package: the connection threshold which is the
#'   minimum signal strength/quality that can be used by a mobile device to connect to an antenna and the probability of
#'   having a two mobile devices
#'
#' @param simFileName the namr of the file used to define a simulation scenario. It is the file that was provided as an
#'   input for the simulation software.
#'
#' @return A list with all the parameters read from the file: start_time, end_time, time_increment, time_stay,
#'   interval_between_stays, prob_sec_mobile_phone, conn_threshold.
#'   

#' @import xml2
#' @export
readSimulationParams <- function(simFileName) {

  if (!file.exists(simFileName))
    stop(paste0(simFileName, " does not exists!"))
  
  simulation.xml  <- as_list(read_xml(simFileName))
  simulation.xml <- simulation.xml$simulation
  simulation.xml$start_time <- as.numeric(simulation.xml$start_time)
  simulation.xml$end_time <- as.numeric(simulation.xml$end_time)
  simulation.xml$time_increment <-
    as.numeric(simulation.xml$time_increment)
  simulation.xml$time_stay <- as.numeric(simulation.xml$time_stay)
  simulation.xml$interval_between_stays <-
    as.numeric(simulation.xml$interval_between_stays)
  simulation.xml$prob_sec_mobile_phone <-
    as.numeric(simulation.xml$prob_sec_mobile_phone)
  simulation.xml$conn_threshold <-
    as.numeric(simulation.xml$conn_threshold)
  
  return (simulation.xml)
  
}