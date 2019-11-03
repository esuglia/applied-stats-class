# Script Class 11
# Fitting a simple Bayesian model using JAGS
# PLS 298 -- F2019
# Andrew Latimer

install.packages("rjags")
install.packages("rstan")
install.packages("MCMCglmm")
install.packages("brms")

library(rjags)
library(MCMCpack)
library(lattice)
library(dplyr)
library(MCMCglmm)
library(brms)

####################################################
#### Part 1: Refitting the model from last week in JAGS ####

# This was the normal model with known variance for weasel body mass. 
# This is simpler than any model you will ever use, but let's use it for direct comparison, then make it slightly more complex.

w <- read.table("Weasel_data_phenotypes.txt", header=T)
w05 <- filter(w, Year==2005 & Gender=="f") # females this time
sigma <- sqrt(16) # we'e saying we already know the variance and it's 16
n.weasels <- length(w05$Body_mass)

# Recall that the model is:
# y~Normal(mu, sigma^2) where sigma is known, and there is the following prior on mu:
# mu~Normal(42.25, 25) -- prior information that the mean is about 42.25 and 
# our uncertainty is symmetrical and equivalent to an sd of 5.

#### A. Coding the model ####

# Now we can code this model in JAGS.
# The main difference between this and R code is that in JAGS (also BUGS) a loop represents taking the product of the individual probability statements in it. So a loop over the data values means that the likelihood is the product of the individual probability calculations for each data point. So for example, the jags code:
tau <- 1/sigma^2
for (i in 1:n.weasels) {
  y[i] ~ dnorm(mu, tau)
}

# Where y is data and mu and sigma are parameters, tells JAGS (or BUGS) that this model has a normal likelihood that is the product of the probabilities of each of the individual y[i] values. 

# The priors on parameters we are estimating are also stated as probability statements. For this example, we have only one parameter to estimate, and thus only one that needs a prior distribution. For JAGS/BUGS, we write the prior for mu as follows:

mu ~ dnorm(42.25, 0.04)

# NOTE: the reason the second parameter is 0.05 (not 25) is that in BUGS and JAGS, the normal distribution is parameterized with the precision (1/variance) instead of the variance or standard deviation.

# Now we can put the model together. The only other logistical detail is that we have to supply the whole model to JAGS as a text file that starts with the word "model" with the whole model enclosed in curly brackets. 
m <- 
"model{ 
  # Likelihood
  for (i in 1:n.weasels) {
    y[i] ~ dnorm(mu, tau)
  }
  
  # Priors
  mu ~ dnorm(42.25, 0.04)  
    # Note the precision tau is known in this model, so there's no prior for it. 
    # Instead, we can just give it a fixed value (1/variance = 1/16):
  tau <- 1/16
}"

# That's it! 
# To use this model, paste it into its own text file and save it. 
# Or write it to a text file from R:
sink(file="weaselknownvariance.jags.txt")
cat(m)
sink()


#### B. Running the model ####

  ####   1. Data = the things we already know ####

# To pass the data to JAGS, we need to make it into a list. 
# The elements of the list are values in the model that we already know.
# These elements of the list also have to have the exact names as in the model code. 

# Q: For this model what are the values that need to be supplied? 
#     And what are they called in the model file? 

jags.data <- list(n.weasels=n.weasels, y=w05$Body_mass)
jags.data

  ####   2.  Initial values = all the other parameters -- what we don't know ####

# We also need to supply initial values for all the model parameters
# which we don't know and which the model will estimate. 
# This is like supplying initial values to an optimizer for finding 
# maximum likelihood values. 

# When running multiple chains, we need to supply a list of initial values for each. For this example, let's use 3 chains, so we need a list of three lists of initial values. 

jags.inits <- list(list(mu=10), list(mu=50), list(mu=150))

  #### 3.  Compiling the model ####

# We tell JAGS what file the model is in and give it the data and initial values. We can also tell it how many chains to create and run, and how long to allow the proposal distribution to adapt. This is the equivalent of setting the "tune" parameter in last week's script, but here JAGS tunes the proposal distribution itself for us. 

weasel.jags <- jags.model(file="weaselknownvariance.jags.txt", data=jags.data, inits=jags.inits, n.chains=3, n.adapt=1000)

# Here is where you almost always will get at least one error message. Most commonly, this happens because of a typo or other syntax error in the model file, or because there is a mismatch between the list of data and/or initial values and the model file. 

#### 4.  Fitting or "updating" the model  ####

update(weasel.jags, n.iter=2000) # run the chains for some number of iterations to allow it to adapt and burn in

# Then update again, saving the samples, to check how it looks
weasel.samples <- coda.samples(weasel.jags, c("mu"), n.iter=1000)
plot(weasel.samples, col=c("black", "blue", "gold"))
acfplot(weasel.samples) # autocorrelation plot for each chain
gelman.diag(weasel.samples) # diagnostic of convergence; the best score is 1, and we can live with <1.2
raftery.diag(weasel.samples)
geweke.diag(weasel.samples)
summary(weasel.samples)


# Now generate a bunch of samples for final parameter estimation.
# If you decide to thin the chain, that can be done now. Here it doesn't seem necessary at all, but we can do it just to see how it works. 

weasel.samples <- coda.samples(weasel.jags, c("mu"), n.iter=10000, thin=1) # now it will save only every 5th sample from the chains and throw the rest away
plot(weasel.samples, col=c("black", "blue", "gold"))
densityplot(weasel.samples)
summary(weasel.samples)
HPDinterval(weasel.samples)

# For comparison, results from a closed-form Gibbs sampler that samples directly from the analytically derived posterior for the mean of a normal distribution with known variance:
w.gibbs <- MCnormalnormal(w05$Body_mass, sigma2=16, mu0=42.25, tau20=25, mc=2000)
plot(w.gibbs)
summary(w.gibbs)


#############################################
#### Part 2 -- Simple linear regression  ####

# You may be thinking that was a lot of work to do the same thing we could do with one line of code using MCMCpack. Yes. But the nice thing is that now with the same infrastructure we can make more complex models. 

# This one's still simple, but starting to be interesting: a linear regression, this time with little or no prior information. 

tree <- read.table("treedata_CO2.txt", header=T)

# There are some missing values in one of the columns of the data set -- for simplicity let's remove those rows. 
tree <- tree[complete.cases(tree),]
dim(tree)

# These data are from the lab manual that accompanies J.S. Clark (2007), "Statistical computation for environmental sciences in R". 
# The data are for growth of trees in a pine plantation near Duke University in which some of the trees were exposed to elevated levels of CO2 for a decade (a FACE experiment). The columns are:
# treeid = tree ID code
# diam.initial = diameter in the first year of the experiment
# diam.final = diameter 10 years later
# treatment = 0/1 where 0 is ambient and 1 is elevated CO2 level

# In this regression we are estimating the effect of elevated CO2 on tree growth. 
# So first let's create a variable for absolute growth over the decade:
tree$growth <- tree$diam.final - tree$diam.initial

# So the model is:
# growth~Normal(mu, tau)
# mu = intercept + slope*treatment 
# plus some priors on the unknown parameters beta (regression coefficients) and tau (precision)

# I typed in a model corresponding to this in the file "treereg.jags".
# Q: BUT before opening this, think about how to write the code for this model. 

  # What will the likelihood include? 
  # Which parameters in the likelihood will need priors?
  # What distributions would be appropriate for those parameters?
  # Note for tau, you can put a prior on the variance or standard deviation if this seems more intuitive, then calculate the precision parameter (tau) from that. 

# Now we go through the same steps as above, but for this more complex model. 

  #### 1.  Organize the data  ####
jags.data <- list(n.trees=dim(tree)[1], treatment=tree$treatment, growth=tree$growth)

  #### 2.  Set initial values for 3 chains ####
# I'm supplying the first ones here, you can fill in the rest. 
jags.inits <- list(list(beta=rnorm(2, 0, 2), sigma2=1), ??? )

  #### 3.  Compile the model and tune the chains ####
tree.jags <- jags.model(file="treereg.jags.txt", data=jags.data, inits=jags.inits, n.chains=3, n.adapt=1000)

  #### 4.  Run the model ####
update(tree.jags, n=1000)
tree.samp <- coda.samples(tree.jags, c("beta", "sigma2"), n.iter=1000, thin=1)
plot(tree.samp)
gelman.diag(tree.samp)
acfplot(tree.samp)
  # a little autocorrelation in the chains for the intercept and slope, nothing that looks like a problem. 

# Run the chains out longer for inference
tree.samp <- coda.samples(tree.jags, c("beta", "sigma2"), n.iter=10000, thin=1)
plot(tree.samp)
gelman.diag(tree.samp)
acfplot(tree.samp)
summary(tree.samp)
  # Here again the Gelman diagnostic suggests the model has converged, and there is little autocorrelation so we don't really need to thin, but we can to it to remove all autocorrelation -- the model runs fast, so why not?

# Run the model again thinning by 5 to retain every 5th value, and saving 2,000 total samples per chain, then recheck these diagnostics. How does it look? 

#### Calculating model DIC ####

# We can also calculate the DIC (deviance information criterion) for this model:
dic1 <- dic.samples(tree.jags, n.iter=5000, thin=5)
dic1 

# For models that have not converged, the DIC value can drift around a lot. To be extra-cautious we can double check by getting the DIC value again and comparing. Note that very small changes in the value are normal, just a consequence of randomness of the MCMC samplers. 
dic.test <- dic.samples(tree.jags, n.iter=5000, thin=5)
dic.test 


# We can compare these results to the classical linear model:
summary(tree.samp)
summary(lm(growth~treatment, tree))
# Q: Any differences? Given the priors used, would you expect any difference? 

# There are R packages that fit particular Bayesian models using MCMC directly in R. A package that fits many mixed models is MCMCglmm. This allows specifying priors, but if you don't specify them, it uses default weak priors on the regression coefficients, random effects, and variance parameters. 
mc1 <- MCMCglmm(growth ~ treatment, random = NULL, data = tree, family = "gaussian")
plot(mc1)
summary(mc1)

# Another, which we will come back to, is the brms library, which interfaces with stan to fit models. 
brms1 <- brm(growth ~ treatment, data = tree, family = "gaussian", set_prior("normal(0,100)"), chains = 2)
summary(brms1)
plot(brms1)
marginal_effects(brms1)
