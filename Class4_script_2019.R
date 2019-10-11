# Script for Class 4
# PLS 298 F2019
# Andrew Latimer


library(bbmle) # Ben Bolker's library of mle functions
library(emdbook) # contains the myxomitosis data set we will use
library(MASS); library(lattice)

#### PART 1: To work through together as a class with projector ####
#   Implementing Bolker's tadpole example of a binomial likelihood

tadpole <- read.table("tadpole.dat")
head(tadpole)
str(tadpole)
# This data set contains results of trials of tadpoles exposed to predators. 
# The 2 relevant columns are: Kill = number killed, and Exposed = number in the trial.

# It's always useful to look at the data.
barplot(tadpole$Kill/tadpole$Exposed, ylim=c(0, 1), xlab="Trials", ylab="Proportion killed")

# Define the negative log-likelihood for the binomial distribution. 
# As in the reading, the binomial likelihood in R is: 
# dbinom(prob=p, k, size=N) Where k is the vector of successes observed,
# N is the number of trials, and p is probability of success, the parameter we're estimating. 

binomNLL1 <- function(p, k, N) { 
  return(-sum(dbinom(prob=p, k, size=N, log=T)))
}
# Q why do you think we want to work on the log scale right away? 
# Alternatively, we could calculate the individual likelihood values and multiply, then log-transform the result.

# We can quickly test to make sure this function is working by supplying a few values for p
binomNLL1(0.5, tadpole$Kill, tadpole$Exposed) # this p is obviously an overestimate since max kills is 5
binomNLL1(0.1, tadpole$Kill, tadpole$Exposed) # this might be closer to correct
# Comparing these two results, does the number change in the direction you expect? What does this indicate?

# The basic R function for optimizing the parameters of a function is optim()
# We can use it to fit the binomial likelihood to data
opt.out <- optimize(f=binomNLL1, interval=c(0,1), N=tadpole$Exposed, k=tadpole$Kill)
# Compare to what we know is the actual ML estimate: 
sum(tadpole$Kill)/sum(tadpole$Exposed)
opt.out

# We could alternatively use Ben Bolker's user-friendlier ML function mle2()
# For the binomial likelihood, it already knows the likelihood function:
mle.out <- mle2(Kill~dbinom(prob=p, size=10), data=tadpole, start=list(p=0.5))
mle.out

# The reason to work with log-likelihoods becomes clear if we look at the actual value of the likelihood here:
logLik(mle.out)
exp(-26.85)
# How would we interpret this probability? Does it matter that it's low? Relative to what? 


#### PART 2: Working in groups  ####
#           Implementing Bolker's myxomitosis example of a gamma likelihood

data(MyxoTiter_sum)

# simplify the data set to just one strain ("grade") of virus
myxdat <- subset(MyxoTiter_sum, grade==1)

# define the Gamma negative log-likelihood
gammaNLL1 <- function(shape, scale) { 
  return(-sum(dgamma(myxdat$titer, shape=shape, scale=scale, log=T)))
}

# We can look up the gamma distribution and see what its moments are to get starting values
shape.start <- mean(myxdat$titer)^2 / var(myxdat$titer)
scale.start <- var(myxdat$titer) / mean(myxdat$titer)

m1 <- mle2(gammaNLL1, start=list(shape=shape.start, scale=scale.start), trace=T)

# Would the normal distribution be a better fit for these data? 
# we could fit the normal and use AIC to compare. 
m2 <- mle2(titer~dnorm(mean=mu, sd=sigma), data=myxdat, start=list(mu=10, sigma=1))

# Or it's more straightforward to use this function from the MASS library:
#m2 = fitdistr(myxdat$titer, densfun="normal") 

AIC(m1, m2) # Which distribution is better according to this criterion?
            # Does the penalty for number of parameters matter here? 


# We might be interested in the shape of the log-likelihood around the MLEs. 
# For a one-dimensional likelihood that's simple. 
# For a 2-dimensional likelihood we need to do a conditional slice through the bivariate surface.
# Ben Bolker wrote some functions that make this easier if we have used his mle2() function to fit the likelihood. 

s1 <- slice(m1,verbose=FALSE) # this takes a "slice" through the negative LL surface, setting the other parameters to their MLE values. 
plot(s1)

# We might also want to see a more general picture of the likelihood "landscape". To get that, we can evaluate the likelihood function across a grid of points, then use a 3d plotting function to display it.

# Set up the values for the x and y axes
x.shape <- seq(20, 80, by=0.5)
x.scale <- seq(0.1, 0.15, by=0.0025)
x.grid <- expand.grid(x.shape, x.scale)
names(x.grid) <- c("shape","scale")
head(x.grid) # check what this data fram looks like
z <- rep(0, dim(x.grid)[1])
for (i in 1:dim(x.grid)[1]) z[i] <- gammaNLL1(x.grid[i,1], x.grid[i,2])
z.mat <- matrix(z, nrow=length(x.shape),byrow=FALSE)
persp(x.shape, x.scale, z.mat, theta=10, phi=10, col="light gray")




#### PART 3: Optional exercise to work through on your own ####
#           Fitting a more complex likelihood

data(ReedfrogFuncresp)

# We can go back to the question of predation on tadpoles, but this time, instead of just fitting a distribution, we can ask whether there is a dependence of predation risk on the density of tadpoles. 

# Bolker demonstrates one way to do this: assume that predators respond to density initially by killing more tadpoles, but that they have a maximum predation rate per unit time, so that the per capita risk to tadpoles declines with increasing density. 

# In this model, there is a "handling time" (parameter h) that each predator needs to process each prey item, so the predators have a maximum predation rate regardless of search time. You could choose different predator response functions but let's follow Bolker and use the Holling Type II (i.e. hyperbolically decreasing risk).

# We're still using the binomial probability to model the data, but now Pr(Kill) is a function of prey density. 

binomNLL2 <- function(params, p, N, k) {
  a <- params[1]
  h <- params[2]
  predprob <-   # fill in the functional response here (answer: see Bolker p.183)
  return(-sum( ??? )) # fill in the likelihood function here (ditto)
}

# Then fit the function
# Because passing the parameters is more complicated in this example, I'll give an example of how you could do it. 
Initial <- as.vector(ReedfrogFuncresp$Initial)
Killed <- as.vector(ReedfrogFuncresp$Killed)
opt2 <- optim(fn=binomNLL2, par= c(a=0.5, h=0.02), N=Initial, k=Killed, control=list(trace=6))



# To look at the model fit we can compare to the data:
a.hat <- opt2$par[1] # MLE for a
h.hat <- opt2$par[2] # MLE for h
x.pred <- seq(1, 101,by=2)
predprob <- a.hat/(1+a.hat*h.hat*x.pred)
predprob.med <- qbinom(0.5, size=x.pred, prob=predprob)
predprob.lcl <- qbinom(0.025, size=x.pred, prob=predprob)
predprob.ucl <- qbinom(0.975, size=x.pred, prob=predprob)

plot(Killed~Initial, ReedfrogFuncresp, pch=16, col="green")
lines(x.pred, predprob.med)
lines(x.pred, predprob.lcl, lty=2)
lines(x.pred, predprob.ucl, lty=2)

# To see how good this model is, what could we compare it to? 
# What are a couple ways you could compare it? 

# for some discussion of the Holling type II function and alternatives, here's one discussion:
# http://www.math.unl.edu/~bdeng1/Teaching/math943/Topics/Stability%20of%20Ecosystems/Skal01.pdf



