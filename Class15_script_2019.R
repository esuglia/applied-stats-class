# Script Class 15
# Fitting a simple population dynamics model to time series data
# PLS 298 -- F2019
# Andrew Latimer

library(MuMIn)
library(ggplot2)
library(lme4)
library(rjags)
library(lattice)

#### Part 1: Data setup and glm model fitting ####

# Here is the population data supplied in McCarthy in jags format
jagsdata <- list(N = matrix(c(32,28,29,39,20,24,22,35,20,36,34,
                             25,25,24,24,28,15,30,35,27,31,20,
                             7,11,8,14,11,5,12,6,10,12,19,
                             12,16,8,8,16,11,10,6,4,10,12), nrow=4, byrow=T), N.years=11, N.pops=4)
                

# If you look at the matrix N you'll see it has one row per population, and each population has possum counts for 11 years (the columns). Since we're modeling the population in any given year as a function of the population in the preceding year, that gives us n.years - 1 = 10 values of N to model per population (we have to start with year 2). 
                
# As a start, we could model this without any explicit population dynamics. That is, we could just model each population value with a linear function of the previous population value plus some noise. 
                
# To make this easier, let's do that kind of simple linear model in R.

# First we need to make a data frame (let's call it "possum") with:
#   1) one column "N" of observed population counts for years 2-11, 
#   2) one column "N.prev" of observed population counts in the previous year,
#   3) one column "pop" as a population index (for the 4 populations),
#   4) one column "year" as a year index (years 2:11)
                
                
# If you're having trouble doing this, you can scroll to the bottom of the script. 

# take a look at the population trajectories
xyplot(N~year|pop,  possum, type="l")
# and the associations between population in current versus previous year
xyplot(N~N.prev|pop, possum)
# possibly more interestingly, the change in population versus population level (i.e. is there density dependence?)
possum$deltaN <- possum$N/possum$N.prev               
xyplot(deltaN~N.prev|pop, possum)
qplot(x=N.prev, y=deltaN, data=possum, colour=factor(possum$pop)) + theme_bw()

# Q: Would you say based on these displays that there is evidence for density dependence in this population time series?

# What then would be a reasonable linear model for these data? 
                
# Given these are count data, lets use a Poisson glm. 
 
m1 <- glm(N~N.prev, data=possum, family="poisson")
# It makes sense that there would be a positive association between last year and this year. We could try to account for density dependence by adding a quadratic term, so that if N.prev was really big, it would start to have a negative effect. 
m2 <- glm(N~poly(N.prev, 2, raw=T), data=possum, family="poisson")

# Now we could start adding variation among populations and/or years with lmer. 

# Try making some lmer models with random intercepts for population and year. 

# Q: which of these score better by AICc? What does this suggest about variation among years and populations?
  

#### Part 2: Bayesian population dynamics model ####

# Now, let's say we want to model the population dynamics more explicitly. That could be useful because then we could have more confidence in simulating those dynamics, to assess extinction risk, to determine whether there are big differences in carrying capacity between sites, and so on.
  
# We know that population in year t is correlated with the previous year's, and also that there is a strong signal of density dependence.  The classic, simple way to represent this kind of population dynamics is a logistic growth model:
#     N(t) = N(t)*r*exp(1-N(t-1)/K)
#   The parameters of this model are:
#       r = maximum population growth rate 
#       K = carrying capacity or equilibrium population density

#   We also know, from our EDA and generalized linear mixed models, that there appear to be meaningful differences among years, so there should be some form of environmental stochasticity built into the model, possibly common across all populations. We could approach this problem purely as a simulation challenge: we could simulate lots of population time series with different values of r and K and try to see which match the dynamics of our data, using some qualitative or quantitative criteria for comparison. 

# We already tried the phenomenological approach of looking for associations between population size and previous population size, and also looking for variation among populations and years. We learned some things from that, but this didn't help us learn about population dynamics parameters. 
  
# Now, let's use the flexibility of Bayesian modeling to do both at the same time. That is, we'll build a population dynamics model into the regression. 
  
# To make sure we have the framework in place, let's first do our simple regression in jags. This is contained in the model file: "possum_linear_jags.txt"

possum1 <- jags.model("possum_linear_jags.txt", data=jagsdata, n.chains=3, n.adapt=3000)                
possum1.samp <- coda.samples(possum1, c("a", "b1", "b2"), n.iter=30000, thin=30)
plot(possum1.samp)
gelman.diag(possum1.samp)
dic1 <- dic.samples(possum1, n.iter=10000); dic1

# Now we can change this by plugging the population process model into the likelihood. See the model file: "possum_logistic_jags.txt"

# Q: What's different in this model from the model possum1? 

possum2 <- jags.model("possum_logistic_jags.txt", data=jagsdata, n.chains=3, n.adapt=3000)                
possum2.samp <- coda.samples(possum2, c("r", "K"), n.iter=10000, thin=10)
plot(possum2.samp)
gelman.diag(possum2.samp)
dic2 <- dic.samples(possum2, n.iter=2000); dic2


                
# Despite having a more biological realistic model, this has higher DIC. 

# What are some potentially important factors we are leaving out? Might there be differences in r and K among populations? 

# Let's look again at the population time series, this time together on a single plot: 
dev.off()
plot(1:11, jagsdata$N[1,], col=1, lwd=3, type="b", ylim=c(0, 50), xlab="year", ylab="Population", cex.axis=1.2, cex.lab=1.5)
for (i in 2:4) lines(1:11, jagsdata$N[i,], col=i, lwd=3, type="b")


#### Part 3: Introducing varying K and environmental stochasticity ####

# Based on what we have seen so far (and on what McCarthy does in the book), we could see whether there's evidence that populations have different carrying capacities. A straightforward way to do this is to compare our current model to a new model that allows K to vary by population. 

# Q: Try modifying the possum2 model yourself to include population-specific K's. Compare using DIC. Is there evidence that different pops have different K's? [again, see below for an implementation of this]

# What about the possibility that there is environmental stochasticity that's common across populations -- i.e. good years and bad years? 

# As an exercise we could try adding year random effects to represent this environmental stochasticity. Does this improve the model? 






#### Hints and solutions ####

# Making possum data frame
N.current <- as.vector(jagsdata$N[,2:11])
N.previous <- as.vector(jagsdata$N[,1:10])
pop.index <- rep(1:4, rep(10,4))
year.index <- rep(2:11, 4)
possum <- data.frame(N=N.current, N.prev=N.previous, pop=pop.index, year=year.index)

# Trying random effects models for possum data using lmer
m1 <- glm(N~N.prev, data=possum, family="poisson")
m2 <- glm(N~poly(N.prev, 2), data=possum, family="poisson")
mer.1 <- lmer(N~poly(N.prev,2) + (1|pop), data=possum, REML=FALSE)
mer.2 <- lmer(N~poly(N.prev,2) + (1|year), data=possum, REML=FALSE)
mer.3 <- lmer(N~poly(N.prev,2) + (1|pop) + (1|year), data=possum, REML=FALSE)
AICc(m1, m2, mer.1, mer.2, mer.3)
BIC(m1, m2, mer.1, mer.2, mer.3)

# Diagnostic plots
par(mfrow=c(2,2)); plot(m2)

# Visualizing a comparison to glm
# Year random effect was preferred, so check whether structure by year in glm residuals.
plot(resid(m2)~possum$year, pch=16, col="navy")
# Hard to say much, with only 4 values per year. 
# Compare fitted vs observed for models without and with random effects. 
par(mfrow=c(2, 1))
plot(fitted(m2)~possum$N, pch=16, col="sienna")
plot(fitted(mer.2)~possum$N, pch=16, col="sienna")

# Model with varying K

possum3 <- jags.model("possum_varyingK_jags.txt", data=jagsdata, n.chains=3, n.adapt=3000)                
possum3.samp <- coda.samples(possum3, c("r", "K"), n.iter=10000, thin=10)
plot(possum3.samp)
summary(possum3.samp)
gelman.diag(possum3.samp)
dic3 <- dic.samples(possum3, n.iter=10000)
dic3

# Better? 

# Model with varying K and also year random effects

possum4 <- jags.model("possum_yearRE_varyingK_jags.txt", data=jagsdata, n.chains=3, n.adapt=3000)                
possum4.samp <- coda.samples(possum4, c("r", "K", "sigma.year", "RE.year"), n.iter=10000, thin=10)
plot(possum4.samp)
summary(possum4.samp)
gelman.diag(possum4.samp)
dic4 <- dic.samples(possum4, n.iter=10000)
dic4

# Do the year random effects improve the model? 