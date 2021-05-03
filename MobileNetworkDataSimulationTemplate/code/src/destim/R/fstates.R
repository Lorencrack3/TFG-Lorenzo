# This function returns the filtered states whenever an
# observation is present and predicted states otherwise.

# P : Transition matrix P(will go to i | is in j)
# E : Event observation matrix P(detection event j | is in i)
# e: Sequence of observation events. Since we are taking
# time increase small, it is expected to have mostly missing
# values. The first value is expected to have an observation.

fstates <- function(P, E, e) {
# Allocate the output matrix
  output <- matrix(0, ncol = length(e), nrow = nrow(P))
# Initial steady state calculation
  edecomp <- eigen(P)
# The steady state might be not unique! Fix this later
  i <- which(abs(edecomp$values - 1) < 1e-6)[1]
  state <- edecomp$vectors[,i]
  state <- state / sum(state)

# Now run the filter. It is a loop, it should be written in
# some compiled language for performance.
  for (i in 1:length(e)) {
    if (!is.na(e[i])) {
      state <- E[,e[i]] * state
      state <- state / sum(state)
    }
    output[,i] <- state
    state <- P %*% matrix(state, ncol = 1)
  }
  return(output)
}
