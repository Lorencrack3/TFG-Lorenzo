## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----LoadLibrary, echo=FALSE--------------------------------------------------
library(destim, warn.conflicts = FALSE)

## ----InitModel, dependson='LoadLibrary'---------------------------------------
model <- HMM(5)

## ----PropertiesInitModel, dependson='InitModel'-------------------------------
nstates(model)
ntransitions(model)
transitions(model)

## -----------------------------------------------------------------------------
constraints(model)

## -----------------------------------------------------------------------------
model <- addtransition(model,c(1,2))
model <- addconstraint(model,c(1,2))
transitions(model)
constraints(model)

## -----------------------------------------------------------------------------
emissions(model)
emissions(model)<-matrix(c(0.3, 0.3, 0.7, 0.9, 0.9,
                           0.7, 0.7, 0.3, 0.1, 0.1),
                         nrow = 5, ncol = 2)
emissions(model)

## -----------------------------------------------------------------------------
model <- HMMrectangle(10,10)
ntransitions(model)
nconstraints(model)

## -----------------------------------------------------------------------------
tws <- matrix(c(3.2, 6.1, 2.2, 5.7, 5.9, 9.3, 5.4,
                4.0, 2.9, 8.6, 6.9, 6.2, 9.7, 1.3),
              nrow = 2, ncol = 7)
S <- function(x) if (x > 5) return(0) else return(20*log(5/x))
emissions(model)<-createEM(c(10,10), tws, S)
dim(emissions(model))

## -----------------------------------------------------------------------------
model <- initparams(model)
all(model$parameters$transitions < 1)
all(model$parameters$transitions > 0)
range(constraints(model) %*% c(model$parameters$transitions, -1))

## -----------------------------------------------------------------------------
ntransitions(model)
model <- minparams(model)
rparams(model)

## -----------------------------------------------------------------------------
obs <- c(1,2,NA,NA,NA,NA,7,7)
logLik(model, obs)
model <- fit(model, obs)
rparams(model)
logLik(model, obs)

## -----------------------------------------------------------------------------
dim(sstates(model, obs))
image(matrix(sstates(model, obs)[,4], ncol = 10))

## -----------------------------------------------------------------------------
dim(scpstates(model, obs))
image(matrix(scpstates(model, obs)[,4], ncol = 100), xlim = c(0,1))

