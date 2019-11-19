# Script Class 12
# Fitting a random effects model using JAGS
# PLS 298 -- F2019
# Andrew Latimer

library(rjags)
library(lme4)

#### Part 1: Fitting a complete-pooling model for the radon data ####

radon <- read.csv("radon.csv")

# make a numeric index for county
counties <- unique(radon$county)
county.index <- match(radon$county, counties)
radon <- cbind(radon, county.index)
head(radon)


#### Part 1 -- Adding random intercepts ####
# Open the file "radon_randomintercept_jags.txt" to look at it.
# For comparison, also open "radon_completepool_jags.txt".

# Q: What is different in the random intercepts model, compared to the complete-pooling model. What do those additional parameters represent?

# complete pooling doesn't have a group index
# random intercepts has a group-level model for the randomly varying intercept
# random intercepts also has a prior on group-level sd and precision
# random intercepts has individual level priors 

# You can also take a look at the "radon_nopool_jags.txt" model for comparison. This model estimates county effects, but without pooling via random intercepts (i.e. they are fixed effects). 


# Data
radon2.data <- list(n=nrow(radon), y=radon$y, x=radon$x, county.index=radon$county.index, n.counties=85)

# for this model we can again let jags generate starting values
radon2.inits <- list(list(sigma.a = 1, sigma.y=2), list(sigma.a = 2, sigma.y=1), list(sigma.a = 5, sigma.y=0.4))
radon2 <- jags.model("radon_randomintercept_jags.txt", data=radon2.data, inits = radon2.inits, n.chains=3, n.adapt=1000)
# run it a bit and check convergence
radon2.samp <- coda.samples(radon2, c("b", "sigma.y", "mu.a", "sigma.a"), n.iter=1000)
summary(radon2.samp)
par(mar=rep(2, 4))
plot(radon2.samp)
acfplot(radon2.samp)
gelman.diag(radon2.samp)

# now get more samples, this time monitoring all parameters of interest
radon2.samp <- coda.samples(radon2, c("a", "b", "sigma.y", "mu.a", "sigma.a", "y.hat"), n.iter=5000, thin=5)

# look at the results
radon2.summ <- summary(radon2.samp) 
radon2.stats <- as.data.frame(radon2.summ$statistics)
head(radon2.stats, 100)
y.hat.rows <- grep("y.hat", rownames(radon2.stats))
y.hat <- radon2.stats$Mean[y.hat.rows] # To get fitted values of the model, look at the stats and pull out the rows corresponding to y.hat
resids <- radon$y - y.hat  
# for convenience, we can then put those values into our data frame and examine them
radon <- cbind(radon, y.hat = y.hat, resid = resids)

# some diagnostic plots
qqnorm(radon$resid)
boxplot(resid~county, data=radon)
# Now that random intercepts for county are included in the model, does it look like the residuals differ strongly by county? 

# Get model DIC
radon2.DIC <- dic.samples(radon2, n.iter=5000, thin=5)
radon2.DIC

# For comparison, let's quickly fit the complete-pooling model that has no county effects. 
radon1.data <- list(n=nrow(radon), y=radon$y, x=radon$x)
# for this simple model we can let jags generate starting values
radon1 <- jags.model("radon_completepool_jags.txt", data=radon1.data, n.chains=3, n.adapt=1000)
radon1.samp <- coda.samples(radon1, c("a", "b", "sigma.y"), n.iter=1000)
radon1.DIC <- dic.samples(radon1, n.iter=1000)



#### Questions for part 1 ####

# A) Is the "complete-pooling" model or the random-intercept model (i.e. "partial pooling") preferred by DIC? How many effective parameters does each model have? How does that compare to the number of fixed and random effects in each model? 

# I think we said partial pooling

# [Bonus: fit the no-pooling "radon_nopool_jags.txt" model. How many effective parameters is that model estimated to have? Is it preferred over the random effects model?]

# B) A couple weeks ago, we fitted a random intercept model to these data using lmer. Re-fit the classical equivalent of this Bayesian model using lmer. What parts of the output from summary() of the lmer() model correspond to the parameters of this model? Do the values for those come out the same?



#### Part 3: Random slope ####

# Add a random slope to the radon2 model. 
# To do this, create a new text file containing the text from "radon_randomintercept_jags.txt", and modify it too add a random slope term. 





#### Also: run a regression model using Stan ####

library(rstan)

stan.data <-  list(N=nrow(radon), K = 2, Y=as.vector(radon$y), X=as.matrix(cbind(rep(1, nrow(radon)), radon$x)))
stanmodelcode <- "
data {
	int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] Y;
}

parameters {
  vector[K] beta;
  real<lower=0> sigma;
}

model {
  vector[N] mu; 
  mu = X * beta;
  beta ~ normal(0, 1000);
  sigma ~ normal(0, 10); 
  Y ~ normal(mu, sigma);
}"

fit.stan = stan(model_code=stanmodelcode, data=stan.data, iter=12000, 
           warmup=2000, thin=10, chains=3)

print(fit.stan)
traceplot(fit.stan)



# We can also quickly fit a random effects model with stan using the brms library, which works great and allows us to make some nice visualizations right off the bat: 
library(brms)
fit.brms <- brm(y ~ x + (1|county), data=radon) # use the default priors
summary(fit.brms)
plot(fit.brms)
plot(marginal_effects(fit.brms))
