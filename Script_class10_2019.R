# Script Class 10
# Writing a function to calculate a posterior probability. 
# Using MCMC to sample from this posterior distribution
# PLS 298 -- F2019
# Andrew Latimer

library(MCMCpack)
library(coda)
library(dplyr)


#### Part 1: Make a function to calculate the posterior distribution of a simple model. ####
# (Normal model with known variance)

# Recall from Bayes theorem what the posterior distribution is: the probability distribution of the parameters, given the prior, the data and the model. 

# For this example, let's visit those Polish weasels again. 
w <- read.table("Weasel_data_phenotypes.txt", header=T)

w05 <- w[w$Year==2005 & w$Gender=="f",] # select one year of data, and females this time

#  For this simple example, we know the true variance -- let's say it's 16. (It is unrealistic that we really would know the variance, but the idea is to keep the example simple.)
sigma <- sqrt(16) # sd=4

# The posterior is proportional to likelihood multiplied times the prior. What is the likelihood and what is a prior for this model?

# There is only one unknown parameter in this model (mu = mean body mass), so the posterior is its distribution. 

# Likelihood

# We've seen this before -- it's y~Normal(mu, sigma^2):

y <- w05$Body_mass # the data

# Here is the likelihood, which is the product of the probabilities of all the individual data points:  
# prod(dnorm(y, mean=mu, sd=sigma)) 
# For computation, it's better to work with the log-likelihood:
# sum(dnorm(y, mean=mu, sd=sigma, log=T))
# This is the relative probability of the data given a value for the parameter mu (mean body mass).

# We can put this into a function:
log.likelihood <- function(mu, sigma, y) { return(sum(dnorm(y, mean=mu, sd=sigma, log=T))) }


# Prior 

# We're not sure what the mean is, but let's imagine that the 2002 data from the weasel data set is in fact a small preliminary study we can use as prior information. 
w02 <- filter(w, Year==2002 & Gender=="f") 
#Our prior estimate of the mean is: 
mean(w02$Body_mass) # prior mean estimate = 42.25

# And that we are pretty highly confident the mean is within 10 g of this value. Assuming we think the prior distribution should be symmetrical, we can use a normal distribution: dnorm(mu, mean=42.25, sd=5). 

# This is the relative probability of the current value of the parameter mu, based on our prior information about mu. Using this prior means we are saying we are 95% sure it's within 10 units of the prior mean.

# We can also put this into a function for convenience:
log.prior.prob <- function(mu) { return(dnorm(mu, mean=42.25, sd=5, log=T)) }


# Posterior: We just multiply the prior and likelihood together, given a value of mu.
mu <- 40 # propose a value for mu, the posterior mean

dnorm(mu, mean=42.25, sd=5) * prod(dnorm(y, mean=mu, sd=sigma))
# Or using the functions we defined, the log posterior probability is:
log.prior.prob(mu) + log.likelihood(mu, sigma, y)

# Using these two pieces, it's the next step is to make a function to return the log posterior probabilty for this model, given a value for mu and data y. Create a function to return the log posterior probability: 

log.post.prob <- function(mu, sigma, y) {
  log.post <- # ???
  return(log.post)
}

# Checking that the function works, and comparing some values for mu (remember, we are assuming we know the value of sigma, so that stays the same).
log.post.prob(mu=30, sigma, y)
log.post.prob(mu=40, sigma, y)
log.post.prob(mu=70, sigma, y)

#### Questions Part 1 ####
# Q: Does it make intuitive sense how these two parts (prior and likelihood) inform the posterior density? What happens to the prior probability when we propose a value for mu that is far from the prior mean? 

#Q: Think about setting up a model with normal likelihood but unknown variance. If you were writing a function to compute the posterior distribution for that model, what would you have to include? 


#### Part 2 -- Use Metropolis sampling to fit this model ####

# The library MCMCpack has a bunch of functions for fitting specific models. It also includes a general-purpose, simple Metropolis sampler, MCMCmetrop1R. We can use that to sample from the posterior density we created in Part 1. 

# First let's run it a short time to look at what it's doing.
w.metrop <- MCMCmetrop1R(log.post.prob, theta.init=70, y=y, burnin=0, mcmc=2000, thin=1, tune=0.4, verbose=1, sigma=sigma)
plot(w.metrop)

# Based on this, what do you think should be the minimum burnin length? Here, the "tuning" parameter controls how big the proposal distribution is. A larger tuning number means the sampler will generate lots of big steps. Big steps are more likely to be rejected, so a higher tuning number means more rejected proposals. Big steps are good for exploring parameter space, but also waste a lot of time sitting at one point before finally accepting a new one. 

# You can see that the function tells you what the acceptance rate is. Generally it's good to have an acceptance rate between 0.25 and 0.5. 

####Question Part 2: Rerun the MCMC sampler, trying different tuning values. What tuning value gives you acceptance rate in the "target" range? How does the chain look when it has that tuning value? ####

# Using what you learned about how this chain behaves, now set the tuning value to a reasonable number, generate a longer MCMC sample, removing some initial samples for burnin, and check it again. This time, check how autocorrelated that samples in the chain are.  

acfplot(w.metrop) # autocorrelation in the chain. 
# You can then thin the chain to remove autocorrelation (changing the "thin" parameter). 

# Rerun the sampler to get at least a few thousand thinned samples. 
# Then we can do some more diagnostics.

acfplot(w.metrop) # check this again after thinning

raftery.diag(w.metrop,  q=0.5, r=0.01, s=0.95) # an estimator for how long you need to run the chain to get a certain level of precision in your estimate of a specified quantile of the distribution.

# Finally, assuming you're happy with the diagnostic results, you can fit the model, incorporating the thinning, tuning, and chain length information you learned. 

# Some ways to look at the output:

plot()
summary()
HPDinterval()
