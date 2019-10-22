# Script for Class 6
# Continuing with random effects models. 
# Random slopes and non-nested random effects. 

# PLS 298 F2019
# Andrew Latimer

library(arm)
library(lme4)
library(lattice)


#### Part 1: Adding random slopes to the model ####

radon = read.table("radon.csv", sep=",", header=T)

# From last time, we had the random-intercepts model:

m3 <- lmer(y~x+(1|county), data=radon)

# What if we thought that there might also be differences in the slope among counties? This would mean that the effect of having a basement on radon levels in a house varies among counties (x=0 means the house has a basement, x=1 means it doesn't).

# The model fitting is easy:

m4 <- lmer(y~x+(1+x|county), data=radon)
display(m4)

# To get an idea what this model is doing, we can look at the overall regression line and compare that to the individual regression lines for counties. 
county.coefs <- coef(m4)$county
plot(0, 0, xlim=c(0, 1),  ylim = c(0, 2), col=0, xlab="x", ylab="y") # set up plot
apply(county.coefs, 1, abline, col=gray(0.7)) # add individual county regression lines, using the apply function to run abline() on every row of the county parameter data set
abline(fixef(m4), lwd=3) # add overall regression line

# Or we could look at the distribution of both group-level coefficients.
par(mfrow=c(2, 1))
hist(coef(m4)$county[,1], main="Intercepts")
hist(coef(m4)$county[,2], main="Slopes")

# A useful kind of plot to make in random effects models is one that displays random effect estimates and uncertainty around them. 

# Defining a function to make a simple version of this kind of plot.
#   coef.mean = vector of estimates of the random effects
#   coef.se = vector of estimates of their standard error
#   grand.mean = overall model mean estimate of the parameter
plot.coef.se <- function(coef.mean, coef.se, grand.mean) {
  n.coef <- length(coef.mean)
  plus.se <- coef.mean+coef.se
  minus.se <- coef.mean-coef.se
  plot(1:n.coef, coef.mean, ylim=c(min(minus.se), max(plus.se)), xlab = "Coefficient", ylab="Parameter value", pch=16)
  abline(grand.mean, 0)
  for (i in 1:n.coef) lines(rep(i,2), c(minus.se[i], plus.se[i]), lwd=2, col=gray(0.7))
  return()
}

# Extract the information from the model
coef.mean = coef(m4)$county
coef.se = se.ranef(m4)$county
grand.mean = fixef(m4)
# Plot the random intercepts
plot.coef.se(coef.mean[,1], coef.se[,1], grand.mean[1])
# and slopes
plot.coef.se(coef.mean[,2], coef.se[,2], grand.mean[2])

# Q: What can you observe from these plots about the random intercepts and slopes and their variability? 

# We can display the same information while relating it to sample size. This will show: 1) how sample estimates tend to stabilize at larger n; and 2) how uncertainty about the coefficients also declines with larger n. 
n.county = tabulate(radon$county) # number of samples per county. 

# modified plotting function
plot.coef.vs.sample.size <- function(coef.mean, coef.se, sample.size, grand.mean) {
  n.coef <- length(coef.mean)
  plus.se <- coef.mean+coef.se
  minus.se <- coef.mean-coef.se
  plot(sample.size, coef.mean, ylim=c(min(minus.se), max(plus.se)), xlab = "Sample size", ylab="Parameter value", pch=16)
  abline(grand.mean, 0)
  for (i in 1:n.coef) lines(rep(sample.size[i],2), c(minus.se[i], plus.se[i]), lwd=2, col=gray(0.7))
  return()
}

# Random intercepts again
plot.coef.vs.sample.size(coef.mean[,1], coef.se[,1], n.county, grand.mean[1])
# Random slopes
plot.coef.vs.sample.size(coef.mean[,2], coef.se[,2], n.county, grand.mean[2])

#Many counties have few or no houses without first-floor measurements, so there's little to no information about the county-specific slope. 

colours <- ifelse(tapply(radon$x, radon$county, max)==0, "red", "darkgray")
plot(coef.mean[,2], col=colours, pch=16)  # counties in which no houses without basements were measured are displayed in red
abline(grand.mean[2], 0)
 
# Q How is the model able to estimate a slope for those counties at all?


#### Part 2: Non-nested random effects ####

flow <- read.csv("streamflowdata.csv")

# These data contain a subset of information on stream flow and characteristics at a number of sites in NE California. We'll look at the response variable "OctFlow", October flow rate. The data were collected at the same 19 sites ("site") in 3 successive years ("year"). 

# For analyzing these data, you might want to consider year as nested within site, but for this exercise, the goal is to work with a model that has two sets of random effects that are not necessarily nested or interacting, so we'll model them as separed (or "crossed") random effects.

# Make sure year, a numeric value, is recognized as a factor:
flow$year = as.factor(flow$year)
levels(flow$year)

# Look at the data
boxplot(OctFlow~Stream, flow)
boxplot(OctFlow~year, flow)
xyplot(OctFlow~year|Stream, flow, type=c("p"))

# Fit a random effects model with random intercept for stream
flow.m1 = lmer(OctFlow~(1|Stream), flow, REML=FALSE)
display(flow.m1)

# Now also add a random intercept for year, since there appeared to be a lot of inter-year variation. 
flow.m2 = lmer(OctFlow~(1|Stream)+(1|year), flow, REML=FALSE)
display(flow.m2)

# Q How much of the variation is associated with site, year, and neither one? 

# Compare models
AIC(flow.m1, flow.m2)


### Part 3 -- Simulating grouped data ####

# This is optional, and explores how estimates of group-level variance behave with different sample sizes and different "true" variance parameters.

# Simulate data from a normally distributed process with n individuals divided among k groups
# with individual and group-level error
#   n = sample size in each group (assumes same sample size in all groups)
#   k = number of groups
#   sigma_y = variance of individuals within groups
#   sigma_alpha = variance among groups

sim.grouped.data <- function(n, k, sigma_y, sigma_alpha) {
  groupRE <- rnorm(k, 0, sigma_alpha)
  groupRE.sim <- rep(groupRE, rep(n, k))
  group <- rep(1:k, rep(n, k))
  y.sim <- rnorm(n*k, 2+groupRE.sim, sigma_y)
  sim.data <- data.frame(y.sim = y.sim, groupRE.sim=groupRE.sim, group=group)
  return(sim.data)
}

#### Small among-group variance
d1 <- sim.grouped.data(n=10, k=50, sigma_y=1, sigma_alpha=0.1)
par(mfrow = c(1, 2), pty="s")
plot(d1$y.sim, xlab="data point", ylab="value", cex.lab=1.5, col="goldenrod3", pch=16)
boxplot(y.sim~group, d1, xlab="group", cex.lab=1.5, col="cadetblue")

#### Large among-group variance
d2 = sim.grouped.data(n=10, k=50, sigma_y=1, sigma_alpha=2)
par(mfrow = c(1, 2), pty="s")
plot(d2$y.sim, xlab="data point", ylab="value", cex.lab=1.5, col="goldenrod3", pch=16)
boxplot(y.sim~group, d2, xlab="group", cex.lab=1.5, col="cadetblue")


#### 3. How many groups do we need to estimate an among-group variance? #### (See G&H p.275)
# We can also use this kind of simulation to look at how number of groups
# affects our ability to recover the "true" parameters. 

# First we can fit models to the simulated data with 50 groups
lmer(y.sim~(1|group), d1)
lmer(y.sim~(1|group), d2)

# What happens if we drop the number of groups to 10? 
d3 = sim.grouped.data(n=10, k=10, sigma_y=1, sigma_alpha=0.1)
d4 = sim.grouped.data(n=10, k=10, sigma_y=1, sigma_alpha=2)
lmer(y.sim~(1|group), d3)
lmer(y.sim~(1|group), d4)

# All the way down to 5?  
d5 = sim.grouped.data(n=10, k=5, sigma_y=1, sigma_alpha=0.1)
d6 = sim.grouped.data(n=10, k=5, sigma_y=1, sigma_alpha=2)
lmer(y.sim~(1|group), d5)
lmer(y.sim~(1|group), d6)
# Q: does the amount of group-level variance appear to affect whether it can be estimated with a small number of groups? 

# Q: if you rerun these simulations a few times, how much do the results vary? 
