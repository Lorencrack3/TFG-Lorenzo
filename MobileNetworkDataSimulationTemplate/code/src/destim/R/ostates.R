# This function returns the observed states whenever an
# observation is present.

# P : Transition matrix P(will go to i | is in j)
# E : Event observation matrix P(detection event j | is in i)
# e: Sequence of observation events. Since we are taking
# time increase small, it is expected to have mostly missing
# values. The first value is expected to have an observation.

ostates <- function(P, E, e) {
  # Allocate the output matrix
  output <- matrix(0, ncol = length(e), nrow = nrow(P))
  # Initial steady state calculation
  edecomp <- eigen(P)
  # The steady state might be not unique! Fix this later
  i <- which(abs(edecomp$values - 1) < 1e-6)[1]
  state <- edecomp$vectors[,i]
  state <- state / sum(state)

  # The transition matrix is only used to compute the
  # steady state.
  for (i in 1:length(e)) {
    if (!is.na(e[i])) {
      output[,i] <- E[,e[i]] * state
      output[,i] <- output[,i] / sum(output[,i])
    }
  }
  return(output)
}
