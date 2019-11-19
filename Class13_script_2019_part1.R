# Script Class 13
# GLMMs for count data
# PLS 298 -- F2019
# Andrew Latimer

library(MuMIn)
library(ggplot2)
library(lme4)
library(rjags)
library(MASS)
library(dplyr)


#### Part 1: Foundation -- build a generalized linear mixed model (GLMM) for count data in BUGS/JAGS ####

# Load data extracted from Kery, "Introduction to WinBUGS for Ecologists". The data include counts  ("count1") of hares at 56 sites ("site") in Switzerland over a seventeen-year period ("year"). The data also incude two environmental variables that are thought to be associated with  hare presence/absence and abundance, "elevation" and "landuse". 

#Note that the sites are also grouped into regions ("region"), so it would be possible to build another level of hierarchical structure into the model. But we won't do that in this script. 

h <- read.table("hares.csv", sep=",", header=T)
head(h)

# These are count data, so we might conventionally start with a statistical model based on the Poisson distribution. First, it can be useful to see if the data are overdispersed, or have excess 0's (zero inflation) or both.

# Compare data to a Poisson distribution
qqplot(h$count1, rpois(10000, lambda=mean(h$count1)))

# Compare data to a negative binomial distribution
z <- fitdistr(h$count1, densfun="negative binomial")
qqplot(h$count1, rnegbin(10000, mu=z$estimate[2], theta=z$estimate[1]))

# Q: are the data obviously overdispersed relative to the Poisson distribution? 

# If so, is the overdispersion of counts mainly due to among-site variation? First is there a lot of variation among sites? 
boxplot(count1~site, h)

# Second, if we look site-by-site, are the data at each location approximately Poisson-distributed? 
sites <- unique(h$site)
par(mfrow=c(4,4), mar=rep(2, 4))
for (i in 1:16) { 
  counts = h$count1[h$site==sites[i]]
  qqplot(rpois(10000, lambda=mean(counts)), counts) 
}

# NOTE:
# Actually in this case, much of the among-site variation is due to variation in sampling intensity -- different-sized areas were surveyed in different locations. So we would ideally want to include in the model an offset for area. 

# But for this exercise, let's proceed as if these sites are heterogeneous in some way other than sampling area. This allows us to compare ways of modeling the resulting overdispersion. 

# As a null model, we could just fit a glm with no site effects:
m0 <- glm(count1~1, data=h, family="poisson")

# Next, as we've done before, to add a random intercept for the grouping variable "site". This time, though, we will use a Poisson distribution with a log link function, rather than a Gaussian / normal distribution.

m1 <- glmer(count1~(1|site), data=h, family="poisson")
summary(m1)
qqnorm(resid(m1))
plot(fitted(m1)~h$count1)
plot(resid(m1)~fitted(m1))

# Just including random effects soaks up a lot of the overdispersion in the data. On the other hand, we might do better if we make the model more flexible to deal with further overdispersion. Because we can make any number of variations on these models in JAGS, we will switch over now and compare among some candidate models, including this random effects model. 

# To see how this model looks in BUGS/JAGS, we can turn to the model contained in "hares_poisson_RE_jags.txt"

# To make this model work, we need to make a numerical index variable "site" for JAGS:
site <- match(h$site, unique(h$site))

hares.data <- list(N=dim(h)[1], N.sites=length(unique(site)), count1 = h$count1, site=site)

m2 <- jags.model("hares_poisson_RE_jags.txt", data=hares.data, n.chains=3, n.adapt=1000)

m2.test <- coda.samples(m2, n.iter=500, c("mu.a", "sigma.a"))
plot(m2.test)
acfplot(m2.test)
gelman.diag(m2.test)

# Run again and look at some diagnostics
m2.samp <- coda.samples(m2, n.iter=2000, c("mu.a", "sigma.a", "resid"))
m2.summ <- summary(m2.samp)
m2.summ$statistics
resid.rows <- grep("resid", rownames(m2.summ$statistics))
m2.resid <- m2.summ$statistics[resid.rows,1]
qqnorm(m2.resid)
# Residuals don't have to be exactly normal, but major deviation from normality may indicate that there are data points that are poorly fitted and/or points that have high leverage on the model results. 

# Get posterior samples of fitted values
m2.samp.lambda <- coda.samples(m2, n.iter=2000, c("lambda"))
m2.fitted <- summary(m2.samp.lambda)$statistics[,1]
# Observed vs fitted
plot(m2.fitted~h$count1); abline(0, 1, lwd=2, col="gray")


#### Part 2: Extending the model to deal with overdispersion ####

# In R, if you discover you have excess 0's or overdispersion, there's are currently only limited ways to deal wiht this while also fitting random effects models. 

# In a Bayesin modeling framework, it's relatively straightforward: make the mean of the Poisson distribution into a random variable. This makes the model effectively into a mixture of Poissons. 

# If lambda is modeled as a gamma distribution, then this is, mathematically, equivalent to a negative binomial model. This isn't the only way to approach this -- for example, we can inject overdispersion by modeling variation in lambda with a lognormal distribution. 

# Check out this model: "hares_poisson_gamma_RE_jags.txt"

# Q: What's different here from m2? 

m3 <- jags.model("hares_poisson_gamma_RE_jags.txt", data=hares.data, n.chains=3, n.adapt=1000)

m3.test <- coda.samples(m3, n.iter=500, c("mu.a", "sigma.a", "alpha"))
# Check whether the markov chains have converged, and are mixing well. 

# Then generate more samples to use for model evaluation.
m3.samp <- coda.samples(m3, n.iter=5000, c("mu.a", "sigma.a", "alpha"), thin=5)
plot(m3.samp)
acfplot(m3.samp)
gelman.diag(m3.samp)
summary(m3.samp)


# The poisson-gamma mixture model is the classic, because it has nice analytical properties and is equivalent to the negative binomial model. 

#   But it's just as easy to use alternative ways of increasing dispersion, for example the Poisson-lognormal model, which is easier to interpret because there is a variance parameter controlling level of dispersion. 

# See "hares_poisson_lognormal_jags.txt"

m4 <- jags.model("hares_poisson_lognormal_RE_jags.txt", data=hares.data, n.chains=3, n.adapt=1000)

m4.test <- coda.samples(m4, n.iter=500, c("mu.a", "sigma.a", "sigma.v"))

# Again, evaluate whether the model has converged and generate more samples. 

# Q: are the parameter estimates for the among-plot variance about the same as for the poisson-gamma model? How about the overall mean? 

# Next, compare the models. What do you conclude?

dicsamp2 <- dic.samples(m2, n.iter=1000); dicsamp2
dicsamp3 <- dic.samples(m3, n.iter=1000); dicsamp3
dicsamp4 <- dic.samples(m4, n.iter=1000); dicsamp4


