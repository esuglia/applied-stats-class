# Script Class 9
# PLS 298 -- F2019
# Andrew Latimer

library(ggplot2)
library(dplyr)

#### Part 1: Normal prior and likelihood ####

#  Normal distribution example -- role of precisions (1/variance) as weights

w <- read.table("Weasel_data_phenotypes.txt", header=T)

# This data set contains phenotypic information on weasels captured in Poland in successive years. 
# The data are from here:
# Zub K, Piertney S, Szafrańska PA, Konarzewski M (2012) Data from: Environmental and genetic influences on body mass and resting metabolic rates (RMR) in a natural population of weasel Mustela nivalis. Dryad Digital Repository. doi:10.5061/dryad.54j79463
# Zub K, Piertney S, Szafrańska PA, Konarzewski M (2012) Environmental and genetic influences on body mass and resting metabolic rates (RMR) in a natural population of weasel Mustela nivalis. Molecular Ecology 21(5): 1283-1293. doi:10.1111/j.1365-294X.2011.05436.x

# For this exercise, let's imagine that the goal of this study is to estimate the mean body mass of weasels at this location. For simplicity, let's also pretend that we have only one year of data, 2004, and that we'll look only at the males.

w04m <- filter(w, Year==2004 & Gender=="m")

# We can easily estimate the mean and variance of the distribution of masses
y <- w04m$Body_mass
data.mean <- mean(y)
data.var <- var(y)
n <- length(y) # sample size

# Now, let's say we have prior information about the size of weasels at this site from the previous sampling year, 2003. 
w03m <- filter(w, Year==2003 & Gender=="m")

# We could use that as prior information to inform our estimate. We can use the mean from the study as the prior mean.
prior.mean <- mean(w03m$Body_mass)
# The prior variance represents our uncertainty about that estimate and so has to take into account the actual variance in the data as well as the sample size. The estimate of uncertainty about a population mean is the standard error, so we can use the sample standard error from the previous study. 


prior.sd <- sd(w03m$Body_mass) # ?? [insert here the standard error of the mean for the w03m data]
# prior.sd = 5 #(part b asks you to change the variance of the prior)
# For doing the weighted means of prior and likelihood, we use the variance instead of the standard deviation:
prior.var <- prior.sd^2

# When combining a normal prior with a normal likelihood, the estimate of the mean is a weighted average of the prior and the estimate from the data (McCarthy p25). The weights are based on the variances of the two estimates, and on the sample size. This will be very familiar from our discussion of how partial partial pooling works in random effects. 
post.mean <- (prior.mean/prior.var + data.mean*n/data.var) / (1/prior.var + n/data.var)

# Note that the weights here are based on the inverse of the variance, or the "precision"
# If we convert the variances to precisions, the weighting becomes more obvious:
prior.prec <- 1/prior.var
data.prec <- 1/data.var
post.mean <- (prior.mean*prior.prec + data.mean*n*data.prec) / (prior.prec + n*data.prec)

# The posterior precision is a weighted average of the precisions, weighted by sample size. 
post.prec <- 1/prior.var + n/data.var
post.var <- 1/post.prec

# Putting these together we can show how the prior combines with the data to give a posterior distribution.
qplot(y, geom="blank") + theme_bw() + scale_x_continuous("Body mass", limits=c(75,105)) + scale_y_continuous("Density") +
  scale_color_discrete(name="Distribution", h=c(270, 90)) + 
  stat_function(fun=dnorm, args=list(mean=prior.mean, sd=sqrt(prior.var)), aes(color="Prior")) + 
  stat_function(fun=dnorm, args=list(mean=data.mean, sd=sqrt(data.var/n)), aes(color="Data")) +
  stat_function(fun=dnorm, args=list(mean=post.mean, sd=sqrt(post.var)), aes(color="Posterior"))


#### Questions for Part 1: ####

# a) How has including the prior information affected the posterior distribution? (i.e. how is the posterior different from the distribution based on data alone?)
# they are essentially identical, which means that the prior is uninformative. This makes sense because the uncertainty is high - it's barely giving you any information. The certainty/informativeness of the prior is weighted by the variance (and precision, 1/var, and sd, var^2) and the sample size (n).

# b) What would this plot look like if the previous data collection had had a much lower variance with the same sample size? 
# the prior sd is really high in part a, so the prior is scaled really far back. When you reduce the variance, the prior becaomes less uniform (uncertain) and more informative, with greater weight put at the mean value.

# Note these plots are similar to Figures 1.4 and 1.5 in the McCarthy reading. 


#### Part 2: Beta prior and binomial likelihood ####

#  Binomial distribution example -- role of sample size

# We can use a hypothetical example of the probability that a fishing trawl will catch a dogfish in a certain area of the ocean. 

# The data are a collection of observations representing success (1) or failure (0) for each Bernoulli trial
d <- c(0,0,0,1,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0)
n <- length(d) # number of observations
y <- sum(d) # number of "successes" (i.e. a dogfish was caught)

# If we just use these data to estimate the frequency of successes, this is just y/n. 
y/n

# But again, let's imagine we had some prior information. Let's say we had some previous observations in the same area that suggested a higher frequency of p = 0.5. To include this information in the model, we need to encode this information into a prior on p, the binomial distribution parameter which we are trying to estimate. 

# The parameter p has support (must lie in the interval) [0,1]. The beta distribution is a suitable candidate because it also has support of [0,1], and also has the other large advantage that it is "conjugate" to the binomial distribution: if you combine a beta prior distribution with a binomial likelihood, you get a posterior with a beta distribution. 

# Here are some examples of the beta distribution with different parameter values
qplot(c(0,1), geom="blank") + theme_bw() + scale_x_continuous("p", limits=c(0,1)) + scale_y_continuous("Density") +
  scale_color_discrete(name="Distribution", h=c(90,270)) + 
  stat_function(fun=dbeta, args=list(shape1=1, shape2=1), aes(color="flat")) + 
  stat_function(fun=dbeta, args=list(shape1=3, shape2=3), aes(color="moderate")) +
  stat_function(fun=dbeta, args=list(shape1=21, shape2=21), aes(color="strong"))

# The beta is also nice to use as a prior for a binomial likelihood because its parameters, a and b, have direct interpretations in terms of successes and failures of a Bernoulli trial. a is equivent to number of successes plus 1, and b to number of failures plus 1. So if we actually have data on trials and outcomes from some other source, we can just include those numbers directly in the prior. If not, we can control the strength of the prior in terms of how many data points we want it to represent. If we think our prior information should be given the same weight as our new data, we could choose values for a and b that imply the same sample size as the data set (n=24):

a.prior <- 13 # if p = 0.5 and n = 24, p/2 is 12, then 13 = # successes + 1
b.prior <- 13 # 12 failures + 1 = 13

a.prior <- 4 # if p = 0.5 
b.prior <- 4 # 3 failures + 1 = 4

# The posterior distribution is straightforward to calculate in this case:

a.posterior <- a.prior + y # prior successes plus observed successes

b.posterior <- b.prior + n - y # prior failures plus observed failures

# The posterior is Beta(p|a.posterior, b.posterior), or in R, dbeta(x, shape1=a.posterior, shape2=b.posterior)

# To compare these distributions, we can plot the prior, the data-based estimate (max likelihood estimate), and the posterior. 

qplot(c(0,1), geom="blank") + theme_bw() + scale_x_continuous("p", limits=c(0,1)) + scale_y_continuous("Density") +
  scale_color_discrete(name="Distribution", h=c(90, 270)) + 
  stat_function(fun=dbeta, args=list(shape1=a.prior, shape2=b.prior), aes(color="prior")) + 
  stat_function(fun=dbeta, args=list(shape1=y, shape2=n-y), aes(color="data only")) +
  stat_function(fun=dbeta, args=list(shape1=a.posterior, shape2=b.posterior), aes(color="posterior"))


#### Questions part 2 ####

# What if we thought the prior information should count only quarter as much as the data? Then we'd have only 6 observations, half of which should be successes. 

a.prior <- 4 # if p = 0.5, 3/6 are successes + 1 = 4
b.prior <- 4 # 3 failures + 1 = 4

a.posterior <- a.prior + y # prior successes plus observed successes

b.posterior <- b.prior + n - y # prior failures plus observed failures

qplot(c(0,1), geom="blank") + theme_bw() + scale_x_continuous("p", limits=c(0,1)) + scale_y_continuous("Density") +
  scale_color_discrete(name="Distribution", h=c(90, 270)) + 
  stat_function(fun=dbeta, args=list(shape1=a.prior, shape2=b.prior), aes(color="prior")) + 
  stat_function(fun=dbeta, args=list(shape1=y, shape2=n-y), aes(color="data only")) +
  stat_function(fun=dbeta, args=list(shape1=a.posterior, shape2=b.posterior), aes(color="posterior"))

# What does the posterior look like now? 
# here, we have less data (n is lower; 6 rather than 24); therefore, the prior is less informative than in the previous example.
