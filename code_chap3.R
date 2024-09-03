# packages
library("deSolve")
library("DEoptim")
library("FME")

# virus dynamics model for acute infections
myModel <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # fill in equations
    dU = - myBeta*U*V
    dI = myBeta*U*V - myDelta*I
    dV = myP*I - myC*V
    
    # outputs of the model
    list(c(dU, dI, dV))
  })
}

myStates <- c(U = 10^6, I = 0, V = 10)
myParams <- c("myBeta", "myDelta", "myP", "myC")
myV <- "V"

modelTime <- seq(from = 0, to = 10, by = 0.01)

# viral load titers (in log base 10)
myData <- read.csv("myData.csv")

# cost function
myCostFn <- function(x) {
  parms <- x[1:length(myParams)]
  names(parms) <- myParams
  ymodel <- ode(myStates, modelTime, myModel, 10^parms)
  yMatch <- ymodel[as.character(ymodel[,1]) %in% as.character(myData$time), ]
  nm <- rle(myData$time)$lengths
  x <- myData[, myV] - rep(log10(yMatch[, myV]), times = nm)
  rmse <- sqrt(mean(x^2))
  return(rmse)
}

## Running the optimizing toolbox
# lower and upper bound
lower = log10(c(1e-7, 1e-2, 1e+0, 1e-1))
upper = log10(c(1e-3, 1e+2, 1e+2, 1e+2))

# optimizer setting
myOptions <- DEoptim.control(itermax = 10000, steptol = 50, reltol = 1e-8)

# fit the model by calling the optimizer and cost function
fit <- do.call("DEoptim", list(myCostFn, lower, upper, myOptions))

# visualize the results
bestPar <- fit$optim$bestmem # parameter estimates in log 10 scale
names(bestPar) <- myParams
out <- ode(myStates, modelTime, myModel, 10^bestPar)

png("./plots/fitted_model.png")
plot(out[, "time"], log10(out[, "V"]), type = "l") # in log 10
points(myData$time, myData$V)
dev.off()
###################################################
### Sensitivity functions 
###################################################

fluCost <- function (bestPar) {
  out <- ode(myStates, modelTime, myModel, 10^bestPar)
  return(cost <- modCost(model = out, obs = myData))
}

fluCost(bestPar)$model

Sfun <- sensFun(fluCost, bestPar)
summary(Sfun)

pairs(Sfun,which=c("V"),col=c("blue"))

# plot
plot(Sfun, which = c("V"), type = "l", lwd = 2)

## the sensitivity parameters 
parRanges <- data.frame(min = log10(0.8*(10^bestPar)), max = log10(1.2*(10^bestPar))) 
rownames(parRanges) <- myParams
parRanges 

tout <- modelTime

## sensitivity to V; equally-spaced parameters ("grid") 
solveFlu <- function(pars) { 
  derivs <- function(t,state,pars) { # returns rate of change
    with (as.list(c(state,pars)), {  
      dU = - myBeta*U*V
      dI = myBeta*U*V - myDelta*I
      dV = myP*I - myC*V
      return(list(c(dU, dI, dV))) 
    }) 
  } 
  state <- c(U = 10^6, I = 0, V = 10)
  tout <- seq(from = 0, to = 10, by = 0.01)
  ## ode solves the model by integration ... 
  return(as.data.frame(ode(y = state, times = tout, func = derivs, 
                           parms = 10^pars))) 
} 

## Plot all variables; plot mean +- sd, min max 
mf <- par(mfrow = c(2, 2)) 
Sens <- summary(sensRange(func = solveFlu, parms = bestPar, dist = "unif", 
                          sensvar = c("V"), parRange = parRanges[4,], num = 10)) 
plot(Sens, log = "y", xlab = "Time (days)", ylab = "log10V", 
     main = "Sensitivity to c", mfrow = NULL) 
par(mf)
title("Viral load sensitivity (varying the respective parameter by 20%)", line = -1, outer = TRUE)

###################################################
### Profile Likelihood
###################################################
myProfile <- function(lower, upper, bestPar) {
  pro.ll <-  NULL
  for (v in 1:length(bestPar)) {
    # Create parameter sequence
    tmpl <- seq(lower[v], bestPar[[v]], length.out = 100)
    tmpl <- tmpl[order(tmpl, decreasing = TRUE)[cumsum(1:13)]]
    tmpr <- seq(bestPar[[v]], upper[v], length.out = 100)
    tmpr <- tmpr[cumsum(1:13)]
    pars <- sort(unique(c(lower[v], tmpl, bestPar[[v]], tmpr, upper[v])))
    ppl  <- NULL
    # Run optimization for each, and record the parameters and RMSE
    for (p in pars) {
      DEargs <- list(myCostFn, replace(lower, v, p), replace(upper, v, p), myOptions) 
      fit <- do.call("DEoptim", DEargs)
      ppl <- c(ppl, fit$optim$bestval)
    }
    pro.ll[[v]] <- cbind(pars, ppl)
  }
  return(pro.ll)
}

outProfiles <- myProfile(lower, upper, bestPar)

# form 26 points
png("./plots/profile_likelihood.png")
par(mfrow = c(2, 2))
sapply(1:4, function(x) plot(outProfiles[[x]], xlab = myParams[x], ylab = 'RMSE'))
dev.off()

## bootstrapping
myBoot <- function(numboot = 1000, numpar = 4) {
  # numpar: number of parameters in the model
  # numboot: number of bootstrap samples
  results <- matrix(NA, numboot, numpar)
  original <- myData
  sampling <- function(x) sample(original$V[original$time==x], length(original$V[original$time==x]), replace = 1)
  for (i in 1:numboot) {
    message("Bootstrapping sample ", i)
    tmp <- sapply(unique(original$time), sampling)
    myData <- cbind(original$time, as.vector(tmp))
    DEarguments <- list(myCostFn, lower, upper, myOptions)
    fit <- do.call("DEoptim", DEarguments)
    results[i, ] <- fit$optim$bestmem
  }
  results <- as.data.frame(results)
  colnames(results) <- myParams
  return(results)
}

bootResults <- myBoot()  

## visualize the results
par(mfrow = c(2,2))
# histogram of boot results
sapply(1:4, function(x) hist(bootResults[, x]) )
# confidence interval of each parameter
apply(bootResults, 2, quantile, probs = c(.025,.975))




