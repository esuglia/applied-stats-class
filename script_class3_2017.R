# Script for Class 3
# PLS 298 F2019
# Andrew Latimer

# if necessary install new libraries
install.packages("mnormt")
library(MASS)
library(mnormt)


#### PART 1: Working with probability distributions in R ####

# Let's start with a data set that has an roughly normal probability distribution
skink <- read.csv("skink_svl.csv")

# This data set has sizes (snout-vent lengths=svl) of skinks by sex (M/F)
# These data allow us to fit a normal distribution. 
qqnorm(skink$svl)
mean.skink <- mean(skink$svl)
sd.skink <- sd(skink$svl)

# Plot the fitted normal distribution
plot(skink$svl, rep(0, length(skink$svl)), ylim=c(0, 0.05), xlim=c(40, 120), ylab="probability density")
curve(dnorm(x, mean=mean.skink, sd=sd.skink), from=30, to=120, n=100, add=T, col="slateblue", lwd=2)
# optionally, add a smoothed density estimator to see how that compares
lines(density(skink$svl), col="orange", lwd=2)
# Does this look different? What might explain the difference? 
# 

# dnorm vs qnorm: https://stats.stackexchange.com/questions/156883/difference-between-density-and-probability

# Q: If this is the true distribution of skink sizes, what's the size of a skink in the 5th percentile? 
qnorm(0.05, mean=mean.skink, sd=sd.skink)
# qnorm() tells you what the value of the variable is at a specified quantile.

# Q: If you catch a skink at random, what's the probability it will be bigger than 100mm?
1-pnorm(100, mean=mean.skink, sd=sd.skink)
# pnorm() gives you the value of the cumulative distribution function (CDF) at the value you specify. 

# Q: Why for this question do we want to know 1 - pnorm()?
# pnorm() looks at the area under the curve to left of a point that you specify; since here we are asking for the probability of the skink size being bigger than the number you specify, you have to get the inverse of pnorm() which is 1-pnorm()

# Q: What is the probability that the skink is between 60 and 80 mm? 
pnorm(80, mean=mean.skink, sd=sd.skink) - pnorm(60, mean=mean.skink, sd=sd.skink)

#### Probabilities of data points ####

# As we'll see in the next class, a "Likelihood" is a function that returns the probability of the data, given a model and a value for the parameters of that model.An extremely simple model is the normal distribution for a data set. 

# Returning to the skink data, we can ask how "likely" or probable an observation is, given the values we have assigned to the parameters. 

# Given that the distribution is normal(mean=81.7, sd=13.9),
#    What is the likelihood of observing a skink of 100mm? 80mm?
#    To help answer this, here's a quick picture of this probabilty distribution: 
curve(dnorm(x, mean=81.7, sd = 13.9), from=40, to=120)
# From this picture, is it more probable that a random skink that you catch will be 80mm or 100mm long? 
# To get a numerical anwser, you can use the function dnorm(), where the mean and sd are the parameters of the fitted distribution. 
# dnorm() is the probability density function (PDF) -- it tells you what is the probability density is for any given value of x, and given some values for the mean and standard deviation.
#dnorm(??, mean = ??, sd = ??)


# CHALLENGE/EXTENSION for Part 1: Binomial data

# Note that the normal is not the only distribution we can use in a statistical model. For many data types it would be inappropriate. 
# One simple example is the Bernoulli distribution, which is a binomial distribution with a single trial. 
# Let's say we suspect the sex ratio in a species is tilted toward females -- 60% female. 

# This means that the distribution that describes the true population is:
# Sex of observed individuals ~ Bernoulli(0.6).
# Or in R: 
# dbinom(n.successes, size=1, p=0.6)
# Assuming independent observations, what is the probability of observing 3 females and 0 males? 


# How likely is it to get that same outcome if the sex ratio is 1/1 (50% female)?




#### PART 2: Comparing some common probability distributions ####

# NOTE Part 2 is more of a reference and review. It gives you an opportunity to look at and play with a few common probability distributions.


#### Continuous distributions ####

# Often observations are said to be "heavy tailed" or leptokurtic. This means there's more probability density in the "tails", the parts of the distribution far from the mean. Some probability distributions with this characteristic are the t and Cauchy distributions. 

# We can compare what those look like by fitting them to the skink data and plotting them. We can use the handy fitdistr() function from the MASS library to do this.

cauchy.params <- fitdistr(skink$svl, "cauchy")
cauchy.params

plot(skink$svl, rep(0, length(skink$svl)), ylim=c(0, 0.05), xlim=c(30, 140), ylab="probability density")
curve(dnorm(x, mean=mean.skink, sd=sd.skink), from=30, to=140, n=100, add=T, col="slateblue", lwd=2)
curve(dcauchy(x, location=cauchy.params$estimate[1], scale=cauchy.params$estimate[2]), from=30, to=140, n=100, add=T, col="darkgreen", lwd=2)

# The normal distribution (or the Cauchy, or t) is not a good one for variance parameters. Variances can't be negative, so to model them, we need a distribution that stays positive. 

# One flexible distribution for positive continuous variables it the gamma:

plot(0,0, xlim=c(0, 20), ylim=c(0, 0.75), col=0)
curve(dgamma(x, shape=0.1, rate=0.1), from=0, to=20, col=1, add=T)
curve(dgamma(x, shape=1, rate=0.5), from=0, to=20,col=2, add=T)
curve(dgamma(x, shape=5, rate=4), from=0, to=20, col=3, add=T)
curve(dgamma(x, shape=5, rate=1), from=0, to=20, col=4, add=T)
curve(dgamma(x, shape=5, rate=0.5), from=0, to=20, col=5, add=T)

# A similarly flexible distribution that takes on values only between 0 and 1 is the beta:

plot(0,0, xlim=c(0, 1), ylim=c(0, 3), col=0)
curve(dbeta(x, shape1=1, shape2=1), from=0, to=1, n=100, col=1, add=T)
curve(dbeta(x, shape1=4, shape2=2), from=0, to=1, n=100, col=1, add=T)
curve(dbeta(x, shape1=2, shape2=4), from=0, to=1, n=100, col=1, add=T)
curve(dbeta(x, shape1=2, shape2=0.7), from=0, to=1, n=100, col=1, add=T)
curve(dbeta(x, shape1=0.7, shape2=2), from=0, to=1, n=100, col=1, add=T)
curve(dbeta(x, shape1=0.1, shape2=0.1), from=0, to=1, n=100, col=1, add=T)

#### Discrete distributions ####

# Another common kind of data is count data -- discrete and positive. 
# The most common distribution to use for count data is the Poisson. 
plot(0:10, dpois(0:10, lambda=3), col="blue", pch=16)
# Note that for the Poisson distribution, there is only one parameter which is both the mean and variance. In practice the distribution of count data often is more dispersed that that. So a common alternative for count data that allows the variance to be higher (or lower) than the mean is the negative binomial. 
points(0:10, dnbinom(0:10, mu=3, size=2), col="red", pch=16)

# Q: What is the probability of observing 4 events under those two distributions? The probability of observing 8? 


#### Multivariate normal distribution ####

# A common distribution in statistical modeling for parameters that may be correlated, such as errors and regression parameters. 
#
# Simulate from a bivariate normal distribution where the two variables
# are positively correlated, rho = 0.6

# Set the parameters
rho <- 0.6 # correlation
mu.a <- 0 # mean of A
mu.b <- 0 # mean of B
sigma.a <- 1 # sd of A
sigma.b <- 2 # sd of B
  # with this information, put together the covariance matrix
sigma <- matrix(c(sigma.a^2, rho*sigma.a*sigma.b, 
    rho*sigma.a*sigma.b, sigma.b^2), nrow=2, byrow=T)
x.mv <- rmnorm(1000, c(mu.a,mu.b), sigma)

# Plot the simulated data points
plot(x.mv, xlab="A", ylab="B", cex.axis=1.4, cex.lab=1.5)


# Plot with shading showing density
smoothScatter(x.mv, bandwidth=0.5, xlab="A", ylab="B", cex.axis=1.4, cex.lab=1.5)
points(x.mv, col="cyan2", pch=20)

# Add some vertical lines to mark places to view conditional distibution of B
x <- -2:2
for (i in 1:5) {lines(c(x[i],x[i]), c(-7,7), col=i, lwd=2)}

# plot distribution of B conditional on different values of A (taking vertical slices through the joint density)
x <- -2:2
plot(rep(0, 15), -7:7, col=0, xlim=c(0, 0.3), ylab="Value of B", xlab="Density of B given A", cex.axis=1.2, cex.lab=1.5)
sample.points <-  seq(-7, 7, length.out=100)
for (i in 1:length(x)) { 
  slice <- cbind(rep(x[i],100), sample.points)
  z <- dmnorm(slice, mean=c(mu.a,mu.b), sigma)
  lines(z,sample.points, col=i, lwd=2)
}
legend("topright", c("A=-2","A=-1","A=0","A=1", "A=2"), col=1:5, lty=rep(1, 5), lwd=rep(3, 5), cex=1.5)

# Just for fun, make a 3d plot of the fitted distribution using the persp() function
x.a <- seq(-4, 4, by=0.1)
x.b <- seq(-8, 8, by=0.2)
x.grid <- expand.grid(x.a, x.b)
z <- dmnorm(x.grid, mean=c(mu.a,mu.b), sigma)
z <- matrix(z, nrow=length(x.a))
persp(x.a, x.b, z, expand=0.5, shade=0.3, phi=50, theta=-10, col="light blue", box=T, xlab="A", ylab="B", zlab="density")


