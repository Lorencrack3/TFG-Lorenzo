---
title: "Introduction to destim package"
author: "Luis Sanguiao"
date: "2020-01-24"
output: html_document
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: references.bib  
---



This vignette contains an introduction to *destim*: its main purpose, syntax and some technical details of its internal behaviour. Some basic knowledge about Hidden Markov Models would be useful to understand how the package works, but not essential to follow this vignette.


## Location of devices
This section contains a brief explanation about the intended use of the package and some remarks on the methodology.

### Introduction
This package purpose is to estimate the spatial distribution of the number of devices, which corresponds to an intermediate step in the general methodology framework described in @WP5Deliverable1.3.

The network events are the result of the constant interaction between a mobile device and a telecommunication network. While the information comprised in these events can be quite complex (and rich), since we are mostly interested in geolocation, a space likelihood that summarizes the event is all we need. The transformation from the network events to the likelihood is quite technical and clearly corresponds to the mobile network operators (MNO), see @WP5Deliverable1.3 for more details.

By space likelihood, we mean the likelihood of the event conditioned on the position of the device. It is quite obvious to identify this likelihood with the emission probabilities of a Hidden Markov model.

### The model
So, we propose a Hidden Markov model to describe the movement of the devices. Those models are quite general and flexible, can be made simple or as complex as wished. Unlike usual, we are not going to estimate the emissions probabilities, that as said we consider known.

As we are mainly interested in estimate their location, and observation events are expected to give mostly information about it,location is a natural choice as state of the model. Since we are considering a tessellation of the map, each possible state would be a tile.

While this could be the simpler approach, it is important to note that more complexity in the state-space might be useful. For example, we might want to represent a car moving north in a highway and a car moving south in the same tile (the tile contains both lanes) by different states, as they are expected to go next to different tiles.

So, if we denote by $n$ the number of tiles, we are going to have not less than $n$ states, and possibly more, so let us say we have $O(n)$ states. In a Hidden Markov model, this means that we have $O(n^2)$ transition probabilities to estimate. Of course, this is not viable, so we are going to fix to zero all transition probabilities to non-contiguous tiles.

Note that in practice, we can do this without losing generality, because given an upper bound for speed $V$ and a lower bound $E$ for the distance between non-contiguous tiles, we can set $\Delta t = \frac{E}{V}$ and the *jump* will no longer be possible.

Now, we have $O(n)$ non zero transitions, which is more affordable, but still very expensive. If we want to reduce this complexity, one option is to classify the states in a certain number of classes, and constrain the transition probabilities to be equal for tiles of the same kind. This only makes sense for periodic tesselations, so it is a strong argument to use periodic tesselations better than other possible choices (Voronoi, BSA, etc.). This is not a limitation of the package though, so it is still possible to estimate models based in non-periodic tilings, but $O(n)$ parameters would have to be estimated. In practice, the package allows any linear constraint between the transition probabilities.

So, we can use constraints to reduce the number of parameters as much as wanted. It is a good idea to keep small the number of (free) parameters: on one hand the likelihood optimization becomes less computationally expensive and on the other hand we get a more parsimonious model.

### Fitting the model
Once we have defined an appropiate model for our problem, the next step is to estimate the (hopefully few) free parameters. As has been already stated, emissions are known, so there are no emission parameters to fit. The initial state is either fixed or set to steady state, so the only parameters to fit in practice are the probabilities of transition.

The method used to estimate the parameters is maximum likelihood, and the forward algorithm computes the (minus) log-likelihood. A constrained optimization is then done. Note that EM algorithm is generally not a good choice for constrained problems, so it is not used in the package.

To get the objective function and the constraints, some previous steps are required to reduce the dimension of the search space. Let $P$ be the column vector of transition probabilities. The linear constraints can be represented as the augmented matrix $(A \vert B)$ so that $AP = B$. After a pivoted QR decomposition is done, we have $R \tilde{P} = Q' B$ where $\tilde{P}$ is a permutation of $P$ and $R$ an upper triangular matrix with non decreasing diagonal elements. Moreover we can express $R$ in blocks as:
$$
R = \left(
\begin{array}{c c }
R_{11} & R_{12} \\
0 & 0
\end{array}
\right)
$$
where $R_{11}$ is a full-rank square upper diagonal matrix. Accordingly we can define blocks for $Q$ and $\tilde{P}$:
\begin{align}
Q & = \left( 
\begin{array}{c c}
Q_1 & Q_2
\end{array}
\right) \\
\tilde{P} & = \left(
\begin{array}{c}
\tilde{P}_1 \\
\tilde{P}_2
\end{array}
\right)
\end{align}

Note that $Q_2' B = 0$, otherwise the constraints can not be fulfilled. The variables in $\tilde{P}_2$ are taken as the free parameters, because being $R_{11}$ full-rank, we have $\tilde{P}_1 = R_{11}^{-1}(Q_1'B - R_{12} \tilde{P}_2)$ so $\tilde{P}_2$ determines $\tilde{P}_1$. These free parameters are transition probabilities and fully determine the transition matrix and thus the likelihood. Moreover, the equality constraints have vanished and now we only have to impose that transition probabilities are between zero and one. Obviously those are linear constraints for $\tilde{P}_2$, so all we have to do is a linear constrained non-linear optimization in the same space. The linear contraints may seem to be a lot, but in practice, a good modeling will make most of the constraints equal, as most of the transition probabilities are going to be equal.

Usually, algorithms for constrained optimization will require an initial value in the interior of the feasible region. To get such initial value, the following algorithm is used:

1. Set transition probabilities to independent uniform $(0,1)$ random variables.
2. Now the constraints do not hold, so the closest point in the constrained space is got through Lagrange multipliers.
3. Now some of the probabilities might be greater than one or smaller than zero. Those are set once again to independent uniforms.
4. Repeat steps 2 and 3 till all transition probabilities are between zero and one.

As already stated, the initial state is set to steady if not fixed. The steady state is calculated as the (first) eigenvector that its eigenvalue is close enough to one. This should be enough because these Markov chains are expected to be both irreducible and aperiodic, otherwise we would have strange movement restrictions.

### The outputs
Once the model has been fit, we can estimate the smooth states by means of the forward-backward algorithm. The smooth states are the mass probability function of the states given the observed data (and the model), thus they kind of summarize all the information available for a given time $t$. So one of the outputs of the package are the smooth states, that can be aggregated to get a space distribution of the number of devices as explained in @WP5Deliverable1.3, section 4.2.

The other main output of the package is the posterior joint mass probability function for two consecutive instants $t, t + \Delta t$ of the states. As it is a posterior probability, it is once again conditioned on all the information available, but it is more dynamic because its time reference is now an interval. A possible analogy would be the smooth positions and speeds of a particle. The former would correspond to position and the later to speed.

Both outputs are needed to estimate the target population.

## Syntax and basic usage
This section explains briefly the main functions of the package.

### Modeling

Obviously, the first step is to create a model. In *destim*, the primary model creator is the function HMM.


```r
model <- HMM(5)
```
When the first parameter is a number, it contains the number of states, so we have five states. Since we have not specified a list of transitions, all states transition to themselves and so only five transitions too.


```r
nstates(model)
#> [1] 5
ntransitions(model)
#> [1] 5
transitions(model)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    2    3    4    5
#> [2,]    1    2    3    4    5
```
The transitions with non zero probability are represented by an integer matrix with two rows, where each column is a transition. The first column is the initial state and the second the final state. Of course, the states are represented by an integer number.

Now, let us look at the constraints.

```r
constraints(model)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    0    0    0    0    1
#> [2,]    0    1    0    0    0    1
#> [3,]    0    0    1    0    0    1
#> [4,]    0    0    0    1    0    1
#> [5,]    0    0    0    0    1    1
```
As we have not specified any constraints, one constraint by state is introduced, the sum of the transition probabilities fixed one initial state has to be one. Otherwise, the transition matrix would not be stochastic. In general, the package adjusts this specific kind of constraints automatically.

The constraints are represented as the augmented matrix of a linear system of equations. The transition probabilities must fulfill the equations, with the same order as shown in transitions function. So the first coefficient in each row is for the transition probability of the transition shown in the first column of the matrix of transitions, and so on.

Both transitions and constraints can be specified as parameters when creating the model. It is also possible to add transitions and constraints later.

```r
model <- addtransition(model,c(1,2))
model <- addconstraint(model,c(1,2))
transitions(model)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    1    2    3    4    5
#> [2,]    2    1    2    3    4    5
constraints(model)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]    1    1    0    0    0    0    1
#> [2,]    0    0    1    0    0    0    1
#> [3,]    0    0    0    1    0    0    1
#> [4,]    0    0    0    0    1    0    1
#> [5,]    0    0    0    0    0    1    1
#> [6,]    1   -1    0    0    0    0    0
```
Now it becomes possible to transition from state 1 to state 2, and in consequence, the first constraint is changed, so we have $p_{11} + p_{12} = 1$ instead $p_{11} = 1$. Moreover, we have added an equality constraint: the first transition probability is equal to the second one. In the matrix of transitions we can see that those first transitions are the transition from one to two and from one to one respectively. The new can be seen in the last row of the new constraints matrix. The constraints to add can also be expressed as the rows of the constraints matrix we want to append.

The emission probabilities are also stored in a matrix, where the element $e_{i j}$ is the probability of the event $j$ to happen, knowing that the device is placed in tile $i$. Of course, the number of possible events is expected to be much smaller than the number of actually observed events. It is possible, though, that the number of columns of the emissions matrix matches the number of actually observed events. This way we do not save memory with the matrix, but it allows us to do the estimations even for a continuous space of observable events.

Note that, in the particular, if any possible event corresponds to a column of the matrix of emissions, each row sums to one. This does not happen in general, as the columns do not need to be exhaustive, and in the continous case they do not even can.

As we have not specified the emission probabilities, they are set to NULL by default.

```r
emissions(model)
#> NULL
emissions(model)<-matrix(c(0.3, 0.3, 0.7, 0.9, 0.9,
                           0.7, 0.7, 0.3, 0.1, 0.1),
                         nrow = 5, ncol = 2)
emissions(model)
#>      [,1] [,2]
#> [1,]  0.3  0.7
#> [2,]  0.3  0.7
#> [3,]  0.7  0.3
#> [4,]  0.9  0.1
#> [5,]  0.9  0.1
```
Emission probabilities are expected to be computed separately, so the model is ready to directly insert the emissions matrix.

Of course, in practice models will have many states and will be created automatically. While the purpose of this package is estimation and not automatic modeling (at least for the moment), some functions have been added to ease the construction of example models.

The function HMMrectangle creates a rectangle model. This model represents a rectangular grid with square tiles, where you can only stay in the same tile or go to a contiguous tile. This means that there are nine non zero transition probabilities by tile. Moreover, horizontal and vertical transitions have the same probability, and diagonal transitions also have the same probability (but different to vertical and horizontal).

As obvious, the rectangle model only have two free parameters to fit, but the number of transitions can be very high even for small rectangles.

```r
model <- HMMrectangle(10,10)
ntransitions(model)
#> [1] 784
nconstraints(model)
#> [1] 782
```
A small rectangle of 10x10 tiles has 784 transitions! Fortunatelly, we only need to fit two parameters. Note that the number of constraints plus the number of free parameters agrees with the number of transitions.

A very simple function to create emissions matrices is also provided. It is called createEM and the observation events are just connections to a specific tower. The input parameters are the dimensions of the rectangle, the location of towers (in grid units) and the distance decay function of the signal strength Martijn Tennekes (2018). In this case, each tile is required to be able to connect to at least one antenna, so no out of coverage tiles are allowed. Note that this is not a requirement for the model, just a limitation of createEM.

```r
tws <- matrix(c(3.2, 6.1, 2.2, 5.7, 5.9, 9.3, 5.4,
                4.0, 2.9, 8.6, 6.9, 6.2, 9.7, 1.3),
              nrow = 2, ncol = 7)
S <- function(x) if (x > 5) return(0) else return(20*log(5/x))
emissions(model)<-createEM(c(10,10), tws, S)
dim(emissions(model))
#> [1] 100   7
```

### Model fitting
Once the model is defined, the next step is to fit its parameters, by constrained maximum likelihood. As already said, the optimizer usually requires an initial guess, so the function initparams obtains an initial random set of parameters.

```r
model <- initparams(model)
all(model$parameters$transitions < 1)
#> [1] TRUE
all(model$parameters$transitions > 0)
#> [1] TRUE
range(constraints(model) %*% c(model$parameters$transitions, -1))
#> [1] -1.110223e-16  0.000000e+00
```
All transition probabilities are between zero and one and the constraints hold, but no observed data is used.

The function minparams reduces the number of parameters of the model to the number of free parameters, as already explained.

```r
ntransitions(model)
#> [1] 784
model <- minparams(model)
rparams(model)
#> [1] 0.30651700 0.02107885
```
So, only two parameters are really needed! It is possible to assign values with rparams, as the optimizer does, but some transition probability might move outside the interval $[0,1]$. The optimization process avoids this problem constraining by linear inequalities.

Now the model is ready to be fitted. Of course, the observed events are needed. Since we have the emissions matrix in the model, the observed events is just an integer vector, that refers to the appropiate column of the emissions matrix.

```r
obs <- c(1,2,NA,NA,NA,NA,7,7)
logLik(model, obs)
#> [1] 9.934475
model <- fit(model, obs)
rparams(model)
#> [1] 4.999942e-01 3.702555e-10
logLik(model, obs)
#> [1] 8.862961
```
Despite the name, logLik returns minus the log-likelihood, so the smaller the better. As in the example, it is possible to introduce missing values for the observations as *NA*.

### The final estimations
Finally, the model is ready to produce some estimations. The main outputs of this package are the smooth states and the smooth consecutive pairwaise states (sometimes called $\xi_{i j}$ in the literature).

The function sstates returns the smooth states as a matrix of *number of states* $\times$ *number of observations* (missing values included) dimensions. So each column represents the space distribution in its corresponding time slot.

```r
dim(sstates(model, obs))
#> [1] 100   8
image(matrix(sstates(model, obs)[,4], ncol = 10))
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

The function scpstates returns the (smooth) joint probability mass function for consecutive states, analogous to the usually denoted $\xi_{i j}$ probabilities in the Baum-Welch algorithm (once convergence is achieved). This time, each column represents the space bi-variant distribution matrix as follows: the probability of the consecutive pair of states $(i,j)$ can be found in the row $400(i-1) + j$.

```r
dim(scpstates(model, obs))
#> [1] 10000     7
image(matrix(scpstates(model, obs)[,4], ncol = 100), xlim = c(0,1))
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

Both each row and each column of the previous image are bidimensional spaces, so it is difficult to visualize: we would need four dimensions instead two! Even so, it is easy to see a diagonal pattern, which is coherent with the transition matrix. The transition matrix only allows transitions to contiguous points, so it is *almost* diagonal. In consequence, so are the $\xi_{i j}$ probabilities.

## Some remarks about computational efficiency
The package is not optimized yet, but some possible bottlenecks have been already noticed. In this section we are going to dicuss them.

### Model construction
As has been already stated, the package is not really ready for model construction, except for the function HMMrectangle, which serves only for toy models. While it is a slow function, its speed can be improved greatly using some compiled language. It is also an embarrassingly parallel problem, so it scales well. On the downside, its complexity increases with the number of states, which can be huge, as it is potentially the area of the country. While a device normally moves around a much smaller area, in practice we need to estimate the position for all the devices, so it is expected to be a very computationally demanding step.

### Model initialization and parametrization
The algorithm used to find an inital value for the optimizer, solves a linear system of equations at each step. While the order of the system grows with the number of transitions and constraints, in practice, most constraints state the equality between two transition probabilities. It is not difficult to get rid of one of the transitions and the constraint in those cases, so in practice the process scales very well for parsimonious models. The function in the package might seem a bit slow, but the elimination of constraints should be written in some compiled language and optimized.

Mostly all said for initialization also goes for parametrization, where the QR decomposition is the *slow* step.

### Forward-Backward algorithm
While it is not implemented yet, the transition matrix can now be easily made sparse, which will improve greatly the computational efficiency of the algorithm. The algorithm is expected to scale well both in memory and speed, and can also be made parallel. As the forward part is used also for likelihood evaluation (the complete algorithm is used for smoothing), a highly efficient forward-backward algorithm is quite important. 

### Likelihood optimization
The purpose of all the previous steps is to make easier the task of the optimizer, so it should not be a problem if everything else is fine. Transition probabilities are required to be between zero and one, what in general means O(n) inequality constraints. A well specified parsimonious model, will often have a lot less constraints, as most of them are duplicated. As a reference, a HMMrectangle(20,20) has 3364 transitions and only 5 constraints are needed in practice.

Duplicated rows are eliminated from the matrix of inequality constraints before calling the optimizer, just in case some of the equalities between probabilities are consequence of the constraints of the model, but are not a specific constraint.

Duplicate rows removal for a matrix of floats is done through a small C++ function, using RCpp package. The function is not very efficient, except when there are few different rows, and it will be replaced by a function that sorts the matrix and then removes duplicates in future versions of the package. 

# References

