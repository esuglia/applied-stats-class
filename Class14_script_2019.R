# Script Class 14
# Applications of hierarchical modeling:
# Abundance model with observation error
# PLS 298 -- F2019
# Andrew Latimer

library(rjags)
library(ggplot2)
library(coda)

### Code is from the online materials for
# Kéry, M. (2010) Introduction to WinBUGS for ecologists
# http://www.mbr-pwrc.usgs.gov/software/kerybook/

# I have annotated the code and made a few modifications. We'll follow the code Kéry uses in the reading (ch. 21) to generate data, display it, then analyze it. 

# Outline:
# Specify biological model: E(Abundance|Environment)
# Relationship of true abundance to vegetative cover
# Add random variation around the mean (Poisson)
# Specify observation model: P(Counted|Present) 
# Relationship of probability of counting to veg cover
# Add random variation (Binomial)
# Simulated count data =
#  Simulated (true) Abundance[i] * Observation prob[i]


#####################################
### Beginning of data generation code

n.site <- 200 # number of sites for which we're simulating lizard counts
vege <- sort(runif(n = n.site, min = -1.5, max =1.5)) # generate the explantory variable, vegetation cover (centered)

# A. Biological process simulation
# Specify biological model: E(Abundance|Environment))
#Relationship of true abundance to vegetative cover

# Set up coefficients for the biological relationship between vegetation and abundance. This set of coefficient values creates a unimodal (hump-shaped) response shape. 

alpha.lam <- 2  			# Intercept 
beta1.lam <- 2				# Linear effect of vegetation
beta2.lam <- -2				# Quadratic effect of vegetation
# Generate values for lambda (mean & variance of Poisson). Note this includes an exponentiation because as is usual in Poisson regression, the regression model is on the log scale.
lam <- exp(alpha.lam + beta1.lam * vege + beta2.lam * (vege^2))

# Confirm what the relationship looks like 
qplot(vege, lam, ylab = "Expected abundance") + theme_bw()

# Add random variation around the mean (Poisson)
# Simulate the "true" unobserved abundances from Poisson distribution, the mean of which is determined by our process model.
# If we want to make sure we all get the same simulated data, we can first set the "seed" of the random number generator. 

set.seed(34567)
N <- rpois(n = n.site, lambda = lam)
table(N)				# Actual distribution of abundances across sites
sum(N > 0) / n.site			# Actual occupancy

# Plot the simulated true abundances vs values of the explanatory variable.
plot(vege, N, main = "", xlab = "Vegetation cover", ylab = "Realized abundance")
curve(exp(alpha.lam + beta1.lam * x + beta2.lam * (x^2)), lwd=2, col="darkgreen", add=T)

# B. Observation process simulation

# Specify observation model: P(Counted|Present) 
# Relationship of probability of counting to veg cover

# Set up coefficients for the detection probability model. These include an intercept, which sets the overall detection probability where vege = 0, and a slope that represents the effect of vegetation cover on probability of detection. 
#   This set of coefficient values creates a negative association between cover and detection. 
alpha.p <- 1				# Intercept
beta.p <- -3				# Strong negative linear effect 
det.prob <- exp(alpha.p + beta.p * vege) / (1 + exp(alpha.p + beta.p * vege)) # this is an inverse logit calculation, to relate the regression values back to the [0, 1] scale of the parameter p.
det.prob
# Plot the relationship between detection probability and vegetation cover. 
plot(vege, det.prob, ylim = c(0,1), main = "", xlab = "", ylab = "Detection probability")
# easier to detect in lower vegetation cover (which is intuitive)
# Plot apparent abundance, taking into account both true abundance and detection probability
expected.count <- N * det.prob
plot(vege, expected.count, main = "", xlab = "Vegetation cover", ylab = 
"Apparent abundance", ylim = c(0, max(N)), las = 1)
curve(exp(alpha.lam + beta1.lam * x + beta2.lam * (x^2)), lwd=2, col="darkgreen", add=T) # Actual mean abundance

# Add random variation (Binomial)
# Simulated count data =
#  Simulated (true) Abundance[i] * Observation prob[i]

# Now simulate the data set for analysis. This consists of observed abundance for T replicate counts at each of n.site sites. 
R <- n.site
T <- 5					# Number of replicate counts at each site
y <- array(dim = c(R, T)) # storage matrix for the simulated data.

for(j in 1:T) { # Note: we don't need a second, nested loop for sites because the true abundances (N) and detection probabilities (det.prob) are supplied as vectors which are n.site long.
    y[,j] <- rbinom(n = n.site, size = N, prob = det.prob)
}
y


sum(apply(y, 1, sum) > 0)		# Apparent occupancy (number of sites that appear to be occupied, based on the survey)
sum(N > 0)				# True occupancy

# Convert the matrix of observed abundances to long or tabular format 
C <- c(y)

# Create columns indexing site, and providing vegetation cover values for each site. 
site = 1:R				# ‘Short’ version of site covariate
site.p <- rep(site, T)			# ‘Long’ version of site covariate
vege.p <- rep(vege, T)			# ‘Long’ version of vegetation covariate
d = cbind(C, site.p, vege.p)		# Check that all went right
head(d)


# End of data simulation
################################


#### Part 1: "Naive" analysis assuming no detection error ####

# One apparently common approach to replicated count data is to take the maximum observed value and assume this comes close to the true abundance at the site. 

max.count <- apply(y, 1, max)
naive.analysis <- glm(max.count ~ vege + I(vege^2), family = poisson)
summary(naive.analysis)

# Q: How well does this capture the coefficient values we actually used to generate the data?

# Calculate predicted values for abundance on the log scale (i.e. not back-transformed to the scale of the observations).
lin.pred <- naive.analysis$coefficients[1] + naive.analysis$coefficients[2] * vege + naive.analysis$coefficients[3] * (vege*vege)

# Display our model results versus the maximum obsereved abundance data, and versus the true mean abundance levels.
par(mfrow = c(1,1))
plot(vege, max.count, main = "", xlab = "Vegetation cover", ylab = "Abundance or count", ylim = c(0, max(N)), las = 1)
lines(vege, lam,col = "darkgreen", lwd = 2)
lines(vege, exp(lin.pred), col = "gold", lwd = 2)
legend("topleft", c("True biological process", "Naive model prediction"), lwd=c(2,2), col=c("darkgreen","gold"))

# the "true biological process" is more accurate b/c it's simulated data - so it predicts that there are actually more lizards at high vegetation than the naive model, and is correcting for detection probability/observation error
# it's easier to detect abundance of lizards at low vegetation cover, so it makes sense that the naive model shifts the mean left (towards higher veg cover) since it's not taking into account detection probability/observation error (and it's based on our data)
# naive model also has a lower abundance at the mean

# Q: Does using the maximum from the repeated counts correct adequately for detection error?

#### Part 2: Analysis using JAGS ####

# Define model
# If you run this code, it will create a text file "BinMix.txt" in your working directory that contains the model code. Kéry runs this in WinBUGS, but we will run it in JAGS. 

sink("BinMix.txt")
cat("
model {

# Priors
 alpha.lam ~ dnorm(0, 0.5) #dunif(-10, 10) # Kéry used uniform truncated priors. Let's try normal priors instead, following Gelman and the Stan team's advice. 
 beta1.lam ~ dnorm(0, 0.5) #dunif(-10, 10) 
 beta2.lam ~ dnorm(0, 0.5) #dunif(-10, 10) 
 alpha.p ~ dnorm(0, 0.5) #dunif(-10, 10) 
 beta.p ~ dnorm(0, 0.5) #dunif(-10, 10) 

# Likelihood
# Biological model for true abundance
 for (i in 1:R) {			# Loop over R sites
   N[i] ~ dpois(lambda[i])
   log(lambda[i]) <- alpha.lam + beta1.lam * vege[i] + beta2.lam * vege2[i]
 }

# Observation model for replicated counts
 for (i in 1:n) {			# Loop over all n observations
   C[i] ~ dbin(p[i], N[site.p[i]])
   logit(p[i]) <- alpha.p + beta.p * vege.p[i]
 }

# Derived quantities
 totalN <- sum(N[])			# Estimate total population size across all sites

}
",fill=TRUE)
sink()

# Put the data together to supply to JAGS
R = dim(y)[1]
n = dim(y)[1] * dim(y)[2]
vege2 = (vege * vege)
sim.data <- list(R = R, vege = vege, vege2 = vege2, n = n, C = C, site.p = site.p, vege.p = vege.p)

# Inits function
Nst <- apply(y, 1, max) + 1 # True counts are random variables to be estimated, so they need starting values
inits <- function(){list(N = Nst, alpha.lam=rnorm(1, 0, 1), beta1.lam=rnorm(1, 0, 1), 
beta2.lam=rnorm(1, 0, 1), alpha.p=rnorm(1, 0, 1), beta.p=rnorm(1, 0, 1))} # function that generates a set of random starting values for one Markov Chain

# Parameters to estimate
params <- c("N", "totalN", "alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta.p")

# MCMC settings
nc <- 3 # number of chains
nb <- 1000 # number of burnin iterations
ni <- 5000 # total number of iterations 
nt <- 5 # how much to thin by

# Compile model and run for some iterations for burnin
BinMix.model <- jags.model("BinMix.txt", data=sim.data, inits, n.chains=nc, n.adapt=1000)
out <- coda.samples(BinMix.model, params, n.iter=ni, burnin=nb, thin=nt)

# Look at a subset of the parameters in the output
params.to.watch <- c("alpha.lam", "alpha.p", "beta.p", "beta1.lam", "beta2.lam", "totalN")
par(mfrow=c(2,3))
traceplot(out[,params.to.watch]) 

# Questions: Do the chains look good? If not, how do you propose to deal with that, and how will you check convergence? Give it a try! 

# Check out the prior specification. What do you notice about the priors on the regression parameters? What happens if you try setting very weak priors on those parameters? 



#### Model evaluation ####

BinMix.summ <- summary(out)
BinMix.summ$stat
# If we look at the bottom of this table we can find the regression parameters. 
# How do these compare to the true parameter values we used to generate the data? 
# To save scrolling, those values are:
# alpha.lam = 2
# beta1.lam = 2
# beta2.lam = -2
# alpha.p = 1
# beta.p = -4

# And how about the empirical patterns? Total N:
sum(N)

# Mean true abundance across all sites
par(mfrow=c(1,1))
plot(N, BinMix.summ$stat[1:n.site], pch=17, col="cyan3", ylab="Estimated abundance")
abline(0,1, lwd=2, col="gray")

# Q: How well did the model estimate the unobserved, true abundances? 


# Let's look at predicted abundance vs maximum observed and true mean abundance at a subset of sites. 

# Pick 16 sites to check on
sel <- sort(sample(1:n.site, size = 16))
sel

# display distribution of latent variable N at these sites and compare to maximum observed counts and to true abundance values 
par(mfrow = c(4,4), mar=rep(2,4))
for (i in 1:16) {
  hist(out[[1]][,sel[i]], col = "grey", main = paste("Site", sel[i]), xlab = "")
  abline(v = Nst[sel[i]]-1, lwd = 3, col = "red") # maximum observed value
  abline(v = N[sel[i]], lwd = 3, col = "blue") # true unobserved value
}
  

# Plot to compare the naive and hierarchical model predictions to the "true" abundance.
par(mfrow = c(1,1))
# extract posterior mean values of regression coefficients
alpha.lam.est <- BinMix.summ$statistics["alpha.lam","Mean"]
beta1.lam.est <- BinMix.summ$statistics["beta1.lam","Mean"]
beta2.lam.est <- BinMix.summ$statistics["beta2.lam","Mean"]
alpha.p.est <- BinMix.summ$statistics["alpha.p","Mean"]
beta.p.est <- BinMix.summ$statistics["beta.p","Mean"]
plot(vege, N, main = "", xlab = "Vegetation cover", ylab = "Abundance", las = 1, ylim = c(0,max(N)))
lines(sort(vege), lam[order(vege)], col = "darkgreen", lwd = 2)
points(vege, exp(lin.pred), type = "l", col = "gold", lwd = 2)
BinMix.pred <- exp(alpha.lam.est + beta1.lam.est * vege + beta2.lam.est * (vege*vege))
lines(vege, BinMix.pred, col = "blue", lwd = 2)
legend("topleft", c("True biological process", "Naive model prediction", "Hierarchical model prediction"), lwd=c(2,2), col=c("darkgreen","gold", "blue"))


#### More questions to think about: ####
# 1) Do you think it would affect how the model would perform if we had fewer sites but more replicates per site? 
# 2) Would the answer depend on how high the true detection rate is? 

