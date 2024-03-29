## Elena Suglia
## Homework 1
## PLS 298, Applied stats modeling

library(tidyverse)
library(dplyr)
library(lattice)
library(vioplot)
library(lme4)
library(arm)
library(bbmle) # Ben Bolker's library of mle functions
library(MASS)

# Question 1 ----

#d = read_table2("CO2_HW1.txt", col_names = TRUE)
#d = as_tibble(d)

d = read.table("CO2_HW1.txt", header = TRUE)
# Check the loaded data
head(d) # looks like it loaded correctly
class(d) # data.frame
str(d) # structure of the dataframe and data types in each column

# Check how the data is distributed
hist(d$logconc)
stem(d$logconc)
qqnorm(d$logconc, main = "Normal Q-Q plot", xlab = "Theoretical quantiles", ylab = "Sample quantiles")
# data does not appear to follow a normal distribution; skewed to the left based on histogram or perhaps overdispersed based on Q-Q plot

ggplot(d, aes(x=Type, y=logconc)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), aes(fill=Type))
boxplot(logconc~Type, d, xlab = "Location")
# these plots tell us that there is virtually no difference in spread of log concentration of CO2 between the 2 sites

ggplot(d, aes(x=Type, y=uptake)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), aes(fill=Type))
boxplot(uptake~Type, d, xlab = "Location")
# however, there is definitely a difference in CO2 uptake rates between sites (higher uptake rates in Quebec)

# strip plots and lattice plots; used in class 1 script but not sure the point and ggplot looks cleaner
stripplot(uptake~logconc|Type, data = d, jitter.data = T, xlab="log CO2 concentration")
stripplot(uptake~logconc, group = Type, data = d, jitter.data = T, xlab="log CO2 concentration")

ggplot(d, aes(x=logconc, y=uptake, group = Type, color = Type)) +
  geom_point() +
  geom_smooth(span=2) +
  theme_classic() +
ggtitle("log CO2 concentration vs CO2 uptake rate") +
  xlab("log CO2 concentration") +
  ylab("CO2 uptake rate")
  #facet_wrap(~Type) #can also add this line look at each site side by side rather than together on the same graph

# higher log CO2 concentration appears positively correlated with uptake rate, and again, there are higher overall uptake rates in Quebec. Here, you are able to see that there appear to be different slopes in responses to log CO2 concentration between the sites: Quebec has a higher response rate; i.e. in Quebec plants increase their CO2 uptake rates more in response to an increase in log CO2 concentration than Mississippi plants do. Other ways to describe it: The relationship between logconc & uptake depends on location. The association between the response and the explanatory variables here depends on the level of a third variable.

coplot(uptake~logconc|Type, d)
# just wanted to practice using coplot function; in reality this is more useful when you have additional explanatory variables that take on more than two values. 

# Need to transform to meet assumption of normality?

# Write down the assumptions of linear regressions here and make sure they are met
# Validity - seems to be a valid test
# Linearity - relationships look pretty linear

## Fit models ----

### Let's compare model fit between several different models: one that has one explanatory variable (uptake rate), one that accounts for location, and one that has an interaction between uptake rate and location.

# Only looks at relationship between log concentration CO2 and CO2 uptake rate

m1 = lm(uptake~logconc, d)
display(m1)

# Accounts for location
m2 = lm(uptake~logconc + Type, d)
display(m2)

# Includes interaction between uptake rate and location
m3 = lm(uptake~logconc*Type, d)
display(m3)

## Assumptions of a linear regression model

# Validity - our question is simply about the effects of one variable on the other; so in this case the regression model seems to be a valid test of this
# Linearity - relationships look fairly linear on the scatterplot
# Independence of errors
# Equal variance of errors
# Normality of errors

# Test for normality and equal variance of errors:

par(mfrow=c(2, 1), mar=rep(3,4), mgp=c(2,1,0))
plot(m3, which=1:2)

# Interpreting the coefficients of the model ----

# Easier to center and scale the explanatory variables before trying to interpret them

d.center <- mutate(d, logconc = scale(logconc, center=TRUE, scale=TRUE))

### Then we can repeat the same regression and compare 

m4 <- lm(uptake~logconc*Type, d.center)
display(m3)
display(m4)

ggplot(d.center, aes(x=logconc, y=uptake, group = Type, color = Type)) +
  geom_point() +
  geom_smooth(span=2) +
  theme_classic() +
  ggtitle("log CO2 concentration vs CO2 uptake rate") +
  xlab("centered log CO2 concentration") +
  ylab("CO2 uptake rate")

#d.center$logconc <- d$logconc - mean(d$logconc, na.rm=T)

ggplot(d.center, aes(x=logconc, y=uptake, group = Type, color = Type)) +
  geom_point() +
  geom_smooth(span=2) +
  theme_classic() +
  ggtitle("log CO2 concentration vs CO2 uptake rate") +
  xlab("centered log CO2 concentration") +
  ylab("CO2 uptake rate")

## Plot fitted vs observed
betas <- coef(m4) # this extracts the coefficients from the linear model
plot(uptake~logconc, data=d.scaled, pch=16)
abline(betas[1], betas[2], col="blue") # regression line for treatment=0
abline(betas[1]+betas[3], betas[2]+betas[4], col="red") # regression line for treatment=1

# Plotting the regression line plus confidence intervals.

## Couldn't figure out how to do it for a dataset with multiple explanatory variables

# set the range of the explanatory variable to use for displaying the predictions
x.pred <- seq(min(d.center$logconc)*1.1, max(d.center$logconc)*1.1, by=0.05) # generates a sequence of values going just outside the range of observed x values

# The function predict() lets us predict the response variable from various levels of the explanatory variable. 
d.pred <- predict(m4, data.frame(logconc=x.pred), se.fit=TRUE, interval="prediction") 
# Make predictions into a table for easier plotting 
erodium_pred_out <- data.frame(x.pred = x.pred, fit = erodium.pred$fit[,1], lwr = erodium.pred$fit[,2], upr = erodium.pred$fit[,3])

# We can also get intervals around that prediction. We can get "predictive" interval, or a "confidence" interval (see below). 
plot(fit ~ x.pred, erodium_pred_out, type="l", ylim=c(-3, 4), ylab="log(Stem Length mm)", xlab="Days to flower")
lines(lwr~x.pred, erodium_pred_out, lty=2)
lines(upr~x.pred, erodium_pred_out, lty=2)
# plot the data on top of the fit
points(log(stem_length)~days_to_flower, d.scaled)

# Are there more points than you'd expect outside the 95% prediction interval? 

# For comparison, here is a plot showing the confidence interval around the regression line:
erodium.pred <- predict(m4, data.frame(days_to_flower=x.pred), se.fit=TRUE, interval="confidence")  
erodium_pred_out <- data.frame(x.pred = x.pred, fit = erodium.pred$fit[,1], lwr = erodium.pred$fit[,2], upr = erodium.pred$fit[,3])
plot(fit ~ x.pred, erodium_pred_out, type="l", ylim=c(-3, 4), ylab="log(Stem Length mm)", xlab="Days to flower")
lines(lwr~x.pred, erodium_pred_out, lty=2)
lines(upr~x.pred, erodium_pred_out, lty=2)
# plot the data on top of the fit
points(log(stem_length)~days_to_flower, d.scaled)

# Interpreting the coefficients of the model

# intercept = When the data is not centered and scaled, the intercept would represent the predicted log CO2 concentration when CO2 uptake rate equals zero, such as in model m3, where the intercept is -4.40. This is actually more interpretable to me than the intercept that comes out of the model on the centered and scaled data, in this case.

# the coefficient for logconc (2.64): an increase in one unit log concentration of CO2 results in an increase in average uptake rate by 2.64 units, all other variables held constant

# coefficient for TypeQuebec (15.94): being in Quebec will increase your uptake rate by 15.94 units, all else held constant

# coefficient for the interaction between logconc and TypeQuebec: (5.95) this represents the difference in slope for log CO2 concentration rate vs uptake rate,  depending on location

# Question 2 ----

# Load the data set “ecdata_HW1.txt”, which includes some growth and flowering time information on some Erodium cicutarium plants from serpentine and non-serpentine environments. The columns are: 
#  sourceSOILTYPE: soil type of source population, 1 = non-serpentine, 2 = serpentine
#  earlylfno: count of leaves early in the plant’s growth
#  totallfno: count of total leaves at end of experiment
#  ffdate: date of first flowering in days after germination

ec = read.table("ecdata_HW1.txt", header = TRUE)
head(ec)
class(ec)
str(ec)

# Fit a normal distribution to the Erodium ffdate data. Also fit a gamma distribution – does this distribution fit the data better or worse than the normal distribution does? Which is “better” by AIC score, or they both about the same? 

# normal distribution ----

mean.ff = mean(ec$ffdate)
sigma.ff = sd(ec$ffdate)

m1ec = mle2(ffdate~dnorm(mean=mean.ff, sd=sigma.ff), data=ec, start=list(mu=10, sigma=1))
#m1ec = fitdistr(ec$ffdate, densfun = "normal")
summary(m1ec)

# gamma distribution ----

# define the Gamma negative log-likelihood

gammaNLL1 <- function(shape, scale) { 
  return(-sum(dgamma(ec$ffdate, shape=shape, scale=scale, log=T)))
}

# We can look up the gamma distribution and see what its moments are to get starting values
shape.start <- mean(ec$ffdate)^2 / var(ec$ffdate)
scale.start <- var(ec$ffdate) / mean(ec$ffdate)

m2ec <- mle2(gammaNLL1, start=list(shape=shape.start, scale=scale.start), trace=T)
summary(m2ec)

AIC(m1ec, m2ec)
# normal distribution fits better

# Calculate the log-likelihood for the normal distribution at the fitted values of the parameters. Show (graphically or in numbers) that the log-likelihood of the data becomes more negative if you shift the mean parameter value away from its maximum-likelihood value.

# calculate the log-likelihood ----

logLik(m1ec) # -2100.618 # higher likelihood
logLik(m2ec) # -2110.891

m3ec = mle2(ffdate~dnorm(mean=40, sd=sigma.ff), data=ec, start=list(mu=10, sigma=1))
logLik(m3ec) # -2440.377

