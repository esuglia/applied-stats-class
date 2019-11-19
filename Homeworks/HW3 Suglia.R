## Elena Suglia
## PLS 298
## Homework 3

library(tidyverse)
library(rjags)
library(brms)
library(MCMCpack)
library(lattice)
library(dplyr)
library(MCMCglmm)
library(rstan)
library(lme4)

# Question 1
# Using the “litterbags.csv” data set, create and run a model in JAGS that corresponds to this one:

# lmer(N_min_rate~Celastrus + (1|Plot), data=litterbags)

# This data set contains nitrogen mineralization rates at 7 different plots where the treatment was whether or not the invasive liana Celastrus orbiculatus was present or not. 

# In your answer, please include the model code and the means and standard deviations of the intercept, slope, group-level (i.e. among-plot) variance and individual-level (within-plot) variance. Briefly report what you did to check model convergence.

d = read_csv("litterbags.csv")

# model to put in a txt file:
#model {
  # Likelihood

#  for (i in 1:n) {
#    y[i] ~ dnorm(y.hat[i], tau.y)
#    y.hat[i] <- a[Plot[i]] + b*x[i] # regression with intercept for each plot
#  }
  
  # group-level model
#  for (j in 1:n.Plots { # 7 plots
#    a[j] ~ dnorm(mu.a, tau.a) # separate intercepts for every plot, no pooling of information
#  }
  
  #priors
#  mu.a ~ dnorm(0, 0.1) # overall mean
#  sigma.a ~ dnorm(0, 0.5)T(0,) # weakly informative prior on the group-level standard deviation
#  tau.a <- pow(sigma.a, -2) # group-level precision
  
  # individual-level priors
#  b ~ dnorm(0, 0.1) #  prior for slope
#  sigma.y ~ dnorm(0, 0.1)T(0,) # weakly informative prior on the individual-level sd
#  tau.y <- pow(sigma.y, -2) # individual-level precision
#}

# data
length(unique(d$Plot)) # 7 plots
d2.data <- list(n=nrow(d), y=d$N_min_rate, x=d$Celastrus, plot.index=d$Plot, n.plots=length(unique(d$Plot)))

# now, let jags generate starting values
d2.inits <- list(list(sigma.a = 1, sigma.y=2), list(sigma.a = 2, sigma.y=1), list(sigma.a = 5, sigma.y=0.4))
d2 <- jags.model("litterbags_jags_hw3.txt", data=d2.data, inits = d2.inits, n.chains=3, n.adapt=1000)
# run it a bit and check convergence
d2.samp <- coda.samples(d2, c("b", "sigma.y", "mu.a", "sigma.a"), n.iter=1000)
summary(d2.samp)
par(mar=rep(2, 4))
plot(d2.samp)
acfplot(d2.samp)
gelman.diag(d2.samp)

# now get more samples, this time monitoring all parameters of interest
d2.samp <- coda.samples(d2, c("a", "b", "sigma.y", "mu.a", "sigma.a", "y.hat"), n.iter=5000, thin=5)

# look at the results
d2.summ <- summary(d2.samp)[1]
d2.stats <- as.data.frame(d2.summ$statistics)
head(d2.stats, 100)
y.hat.rows <- grep("y.hat", rownames(d2.stats))
y.hat <- d2.stats$Mean[y.hat.rows] # To get fitted values of the model, look at the stats and pull out the rows corresponding to y.hat
resids <- d$N_min_rate - y.hat  
# for convenience, we can then put those values into our data frame and examine them
d <- cbind(d, y.hat = y.hat, resid = resids)

# Some summary statistics:
# The means and standard deviations of:

# - slope = b
# - intercept = mu.a
# - group-level (i.e. among-plot) variance = sigma.a
# - individual-level (within-plot) variance = sigma.y

woo = summary(d2.samp)[1]
woo$statistics[,c('Mean', 'SD')]

# some diagnostic plots
qqnorm(d$resid)
boxplot(resid~Plot, data=d)
# Now that random intercepts for county are included in the model, does it look like the residuals differ strongly by county? 

# Get model DIC
d2.DIC <- dic.samples(d2, n.iter=5000, thin=5)
d2.DIC

# Briefly report what you did to check model convergence.

# For comparison, let's quickly fit the complete-pooling model that has no county effects. 
#d1.data <- list(n=nrow(d), y=d$N_min_rate, x=d$Celastrus)
# for this simple model we can let jags generate starting values
#d1 <- jags.model("d_completepool_jags.txt", data=d1.data, n.chains=3, n.adapt=1000)
#d1.samp <- coda.samples(d1, c("a", "b", "sigma.y"), n.iter=1000)
#d1.DIC <- dic.samples(d1, n.iter=1000)


# Question 2
# Fit a hierarchical (multilevel) model in JAGS using the data sets “immunity.csv” and “patient_age.csv”. The model should include “immune.level” as the response variable, “time” as an individual-level predictor, and “age” as a group-level predictor. 

im = read_csv("immunity.csv")
pa = read_csv("patient_age.csv")

# These data are repeated measures on individual patients, so each patient is a “group” in the data set. Each row in the data set “immunity.csv” is one observation on one patient. Each row in the data set “patient_age.csv” is the age of each patient at time the study began – so this is a group level predictor, with one row of data per patient. In other words, here “age” is analogous to county-level bedrock uranium content in the Gelman & Hill radon example. 

# Hints: You can use the column “patient” in “immunity.csv” to index the random intercept for patient. You can then use the column “age” in “patient_age.csv” in the group-level regression that explains some of the variation in the random intercepts. There is a model like this on page 361 of Gelman & Hill.  

# In your answer, include the model you created. Also report the means and standard deviations of:
# - the slopes of the individual-data-point-level and “group”-level (i.e. patient-level) regressions
# - the individual-data-point-level and “group”-level variance parameters

# model to put in a txt file:
model {
  for(i in 1:n){
    y[i] ~ dnorm(y.hat[i], tau.y)
    y.hat[i] <- a[patient.index[i]] + b*x[i]
  }
  b ~ dnorm(0, .0001)
  tau.y <- pow(sigma.y, -2)
  sigma.y ~ dunif(0, 100)
  
  for(j in 1:n.patients){
    a[j] ~ dnorm(a.hat[j], tau.a)
a.hat[j] <- g.0 + g.1*age[j]
  }
  g.0 ~ dnorm(0, .0001)
  g.1 ~ dnorm(0, .0001)
  tau.a <- pow(sigma.a, -2)
  sigma.a ~ dunif(0, 100)
}

# data
im.data <- list(y=im$immune.level, x=im$time, patient.index = im$patient, age=pa$age, n.patients=length(unique(im$patient)), n = nrow(im))

# now, let jags generate starting values
im.inits <- list(list(sigma.a = 1, sigma.y=2), list(sigma.a = 2, sigma.y=1), list(sigma.a = 5, sigma.y=0.4))
im.model <- jags.model("immunity_jags_hw3.txt", data=im.data, inits = im.inits, n.chains=3, n.adapt=1000)

# run it a bit and check convergence
im.samp <- coda.samples(im.model, c("b", "sigma.a", "sigma.y", "tau.a"), n.iter=1000)
summary(im.samp)
par(mar=rep(2, 4))
plot(im.samp)
acfplot(im.samp)
gelman.diag(im.samp)

# now get more samples, this time monitoring all parameters of interest
im.samp.full <- coda.samples(im.model, c("a", "b", "tau.a", "sigma.a", "sigma.y", "y.hat"), n.iter=5000, thin=5)

# look at the results
im.summ <- summary(im.samp.full)[1]
im.stats <- as.data.frame(im.summ$statistics)
head(im.stats, 100)
y.hat.rows <- grep("y.hat", rownames(im.stats))
y.hat <- im.stats$Mean[y.hat.rows] # To get fitted values of the model, look at the stats and pull out the rows corresponding to y.hat
resids <- im$immune.level - y.hat  
# for convenience, we can then put those values into our data frame and examine them
imres <- cbind(im, y.hat = y.hat, resid = resids)

# Some summary statistics:
# The means and standard deviations of:

# - slope = b
# - intercept = a
# - group-level (i.e. among-patient) variance = sigma.a
# - individual-level (within-patient; individual-data-point-level) variance = sigma.y
# - y.hat?

woo1 = summary(im.samp.full)[1]
woo1$statistics[,c('Mean', 'SD')]

# some diagnostic plots
qqnorm(imres$resid)
boxplot(resid~patient, data=imres)
# Now that random intercepts for patient are included in the model, does it look like the residuals differ strongly by patient? 

# Get model DIC
im.DIC <- dic.samples(im.model, n.iter=5000, thin=5)
im.DIC

