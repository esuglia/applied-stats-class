# Script for Class 16
# PLS 298 F2019
# Andrew Latimer

library(ggplot2)
library(rjags)
library(lattice)


#### Part 1. Random walk time series that includes a time trend ####

# In this example we'll use simulated data and see if we can recover the known variance parameters using a state space model. 

# Set Parameter Values
mu <- 1 # trend (effect of time) for the time series
sigma_proc <- 1 # standard deviation for residual process "error"
sigma_obs <- 0.5 # standard deviation for observation error
nt <- 100 # number of time steps
x0 <- 20 # initial state 

#Simulate Random Walk with trend
x = rep(0, nt)
x[1] = x0 + rnorm(1, mean=0, sd=sigma_proc) + mu
for (i in 2:nt) x[i]= x[i-1] + rnorm(1, mean=0, sd=sigma_proc) + mu # process simulation
y=x+rnorm(nt, sd= sigma_obs) # adding observation error

#View Data
plot(1:nt, x, pch=16, col="darkgray")
points(y, pch=16, col="darkgreen")
legend("topleft", c("True state", "Observed"), col=c("darkgray", "darkgreen"), pch=rep(16, 2))

# State space model

time_trend.jags <-
  "model{

# likelihood for observation / data model
for (i in 1:n) { 
y[i] ~ dnorm(x[i], tau.obs)
}
# likelihood for process/ state model 
for (i in 2:n) {
  x[i] ~ dnorm(mu[i], tau.proc)
  mu[i] <- x[i-1] + trend
}

# priors
x[1] ~ dnorm(0, 0.001)
trend ~ dnorm(0, 0.1)
tau.obs <- pow(sigma.obs, -0.5)
tau.proc <- pow(sigma.proc, -0.5)
sigma.obs ~ dnorm(0, 0.5)T(0,)
sigma.proc ~ dnorm(0, 0.5)T(0,)
}"

sink("time_trend.jags")
cat(time_trend.jags)
sink()

# run the model in jags

jagsdata.tt <- list(y = y, n=nt)
timetrend <- jags.model("time_trend.jags", data=jagsdata.tt, n.chains=3, n.adapt = 1000)
update(timetrend, n.iter=10000)

timetrend.samples <- coda.samples(timetrend, n.iter=10000, thin=10, c("sigma.obs", "sigma.proc", "trend"))

par(mar=rep(1,4)); plot(timetrend.samples)
summary(timetrend.samples)

#### Questions Part 1 ####
# How well does the model estimate the trend? 
# Does it appear to be able to separate observation error from residual process error? 
# Do the estimates get worse if you shorten the time series? 


#### Part 2. State space model of population dynamics #####

# Set Parameter Values
r <- 0.3 # initial (low-population) growth rate
K <- 50
sigma_proc <- 1 # standard deviation for process "error"
sigma_obs <- 0.5 # standard deviation for observation error
nt <- 30 # number of time steps
x0 <- 80 # initial state 

#Simulate population time series
x <- rep(0, nt)
x[1] <- x0*exp(r * (1-x0/K)) + rnorm(1, mean=, sd=sigma_proc)
for (i in 2:nt) x[i] <- x[i-1]*exp(r * (1-x[i-1]/K)) + rnorm(1, mean=0, sd=sigma_proc) # process simulation
y <- x+rnorm(nt, sd= sigma_obs) # adding observation error

#View Data
plot(1:nt, x, pch=16, col="darkgray")
points(y, pch=16, col="darkgreen")
legend("topleft", c("True state", "Observed"), col=c("darkgray", "darkgreen"), pch=rep(16, 2))

# State space model

logistic_growth.jags = 
  "model{

# likelihood for observation / data model
for (i in 1:n) { 
y[i] ~ dnorm(x[i], tau.obs)
}
# likelihood for process/ state model 
for (i in 2:n) {
x[i] ~ dnorm(mu[i], tau.proc)
mu[i] <- x[i-1]*exp(r*(1-x[i-1]/K))
}

# priors
x[1] ~ dnorm(0, 0.001)
r ~ dgamma(1.5,1.5)
K ~ dunif(1, 1000)
tau.obs <- pow(sigma.obs, -0.5)
tau.proc <- pow(sigma.proc, -0.5)
sigma.obs ~ dnorm(0, 0.5)T(0,)
sigma.proc ~ dnorm(0, 0.5)T(0,)
}"

sink("logistic_growth.jags")
cat(logistic_growth.jags)
sink()

# run the model in jags

jagsdata.lg <- list(y = y, n=nt)
logisticgrowth <- jags.model("logistic_growth.jags", data=jagsdata.lg, n.chains=3, n.adapt = 1000)
update(logisticgrowth, n.iter = 5000)

logisticgrowth.samples <- coda.samples(logisticgrowth, n.iter=10000, thin=10, c("sigma.obs", "sigma.proc", "r", "K"))

par(mar=rep(1,4)); plot(logisticgrowth.samples)
summary(logisticgrowth.samples)

#### Questions Part 2 ####
# What would you say about the model's ability to recover the true population dynamics parameters? How about the variance components? 

# If the model is not doing well on estimating the variances, what could you potentially do to improve the situation? What information would be useful to have? 

# If you had information about the detection rate of the organism, where in the model would you include that? 


#### Part 3: Missing values ####

# It's common to have missing values from particular location x time combinations, or in the case of a simple time series, missing observations. 

# In a Bayesian model you can handle these missing values by setting them up as random variables to be estimated. For example, let's say there are several missing values in the simulated population time series we have been working with.  

# Set Parameter Values
r <- 0.2 # trend (effect of time) for the time series
sigma_proc <- 1.5 # standard deviation for residual process "error"
sigma_obs <- 0.5 # standard deviation for observation error
nt <- 50 # number of time steps
x0 <- 30 # initial state 

#Simulate population time series
set.seed(12345)
x <- rep(0, nt)
x[1] <- x0*exp(r * (1-x0/K)) + rnorm(1, mean=, sd=sigma_proc)
for (i in 2:nt) x[i] <- x[i-1]*exp(r * (1-x[i-1]/K)) + rnorm(1, mean=0, sd=sigma_proc) # process simulation
y <- x+rnorm(nt, sd= sigma_obs) # adding observation error

#View Data
plot(1:nt, x, pch=16, col="darkgray")
points(y, pch=16, col="darkgreen")
legend("topleft", c("True state", "Observed"), col=c("darkgray", "darkgreen"), pch=rep(16, 2))

# Let's get rid of some of the values, to treat them as missing. 
missing.index <- c(6:10, 15, 25, 36:45)
y[missing.index] <- NA

# When there are NA's in the response variable, JAGS/BUGS will automatically put a vague normal prior on them and estimate them. So we can just use the same model as before. 

# run model
jagsdata.lg <- list(y = y, n=nt)
logisticgrowth <- jags.model("logistic_growth.jags", data=jagsdata.lg, n.chains=3, n.adapt = 1000)
update(logisticgrowth, n.iter=5000)

logistic.samples <- coda.samples(logisticgrowth, n.iter=10000, thin=10, c("sigma.obs", "sigma.proc", "r", "K", "x"))
out <- summary(logistic.samples)$statistics
round(out[c("r", "K", "sigma.obs", "sigma.proc"),], 2)

# Take a look at the estimates for the state variable
x.mean <- out[grep("x", rownames(out)),1]
x.sd <- out[grep("x", rownames(out)),2]
x.lo <- x.mean - 2*x.sd # approximate credible intervals for x
x.hi <- x.mean + 2*x.sd
plot(x.mean, pch=16, cex=1.5, col=ifelse(is.na(y), "orange", "blue"), ylim=c(min(x.lo), max(x.hi)), xlab="Time", ylab="Population density")
for (i in 1:length(x.lo)) lines(rep(i, 2), c(x.lo[i], x.hi[i]), col="darkgray")
points(y)
abline(out["K",1], 0, lty=2)
legend("bottomright", c("state with y observed", "state with y missing", "y"), pch=c(16, 16, 1), col=c("blue", "orange", "black"))

#### Part 3 questions ####
# How the model is estimating the state values during the times where there are no observations? 
# What happens to the uncertainty about the state variable when one or more observation values are missing? Why? 


