# Script for Class 2 
# PLS 298 -- Applied Statistical Modeling -- F2019

#### Setup ####

library(dplyr)
library(lattice)
library(lme4)
library(arm)
library(ggplot2)
library(viridis)

#### PART 1: Multiple regression ####

## Why might we want or need to include multiple variables in a model?

# a) Include another variable to explain more variation (when no interaction between variables but both are explanatory). 
# Example:
# Predicting the length of petals of Iris versicolor from other morphological measurements 

# In this case, excluding the additional relevant variable results in "underfit". 

plot(Petal.Length~Petal.Width, iris, subset=iris$Species=="versicolor")
plot(Petal.Length~Sepal.Length, iris, subset=iris$Species=="versicolor")

# Here's a useful plot similar to xyplot that plots the relationship between 2 data columns conditioned on a third. 
coplot(Petal.Length~Petal.Width|Sepal.Length, iris[iris$Species=="versicolor",])
# the relationship between petal.length and petal.width, conditional on the sepal.length. The horizontal bars on top graph (sepal length) from lowest to highest each correspond to a subplot below from low left to top right.

# Does the model get better when we add another explanatory variable? 
summary(lm(Petal.Length~Petal.Width, iris, subset=iris$Species=="versicolor"))
summary(lm(Petal.Length~Petal.Width + Sepal.Length, iris, subset=iris$Species=="versicolor"))


# b) Include another variable beacuse the association between the response and the explanatory variables depends on level of other variables.

# When people write that an analysis "controlled for" other variables, they often simply mean that other variables were included in the analysis. They didn't control those other variables in the experimental sense. Instead, the claim is that the analysis is able to account for dependence in the relationship between the response and the explanatory variable of interest. 

# Example:
# Predicting ANPP from rainfall from the grassland data from last class. We might expect the relationship between productivity and rainfall to depend on temperature. We can explore this graphically, and include temperature as a covariate in a regression.

ANPP <- read.csv("Grassland_NPP_31_site_summary.csv", na.strings = "-9999")
# Let's drop the extreme high value for now
ANPP <- filter(ANPP, ANPP1 < 1000)

# Graphical display
coplot(ANPP1~Rainfall|Temperature, ANPP)
# Regression with rainfall alone
summary(lm(ANPP1~Rainfall, ANPP))
# Adding temperature as a covariate
summary(lm(ANPP1~Rainfall*Temperature, ANPP))
        
# Does the effect of rainfall seem to be stronger or weaker when we incorporate temperature into the model? What's the interpretation of the regression coefficients in this second model? 
# cooler temperatures require less rain to increase NPP (because NPP starts at low levels), but at higher temperatures more rain is required to increase NPP because there is more evaporative demand

# We can also look at this interaction by splitting the data into "high" and "low" temperature locations.
# Maybe the easiest way is to manually create a categorical variable first, then use it to differentiate high and low temperature sites in the plot
ANPP$Temp_cat <- ifelse(ANPP$Temperature>10, "Hi", "Lo") 
ggplot(ANPP, aes(x=Rainfall, y=ANPP1)) + geom_point(aes(color=Temp_cat)) + theme_bw()

# Note this could also be done using temperature as a continuous variable 
ggplot(ANPP, aes(x=Rainfall, y=ANPP1)) + geom_point(aes(color=Temperature)) + scale_color_viridis() + theme_classic()


#### PART 2: Regression fitting, diagnostics, interpretation ####

# When you have a continuous and a categorical variable, the requirement to sample across the whole range equates to balance.

# As we'll see, models with random effects are more robust to imbalance, but there's no way around needing to have a big enough sample size that is well stratified enough to learn about each level of a factor. 

# (Note here is another reason you might choose to use a random effect -- you can have low sample size in some groups/categories, you just won't learn much about that particular one but you can still generalize across them). 

d <- read.table("Erodium_data_for_class2.csv", sep=",", header=T)

# Data set explanation: 
# These data are from a greenhouse experiment in which widespread invasive plant Erodium cicutarium was planted in high and low-water treatments. This exercise explores the effect of plant phenology (how soon they flower) and watering treatment on plant size (stem length).
# -- treatment (1/0) = low-water treatment (1) or high-water treatment (0)
# -- stem_length = length of longest stem on each plant in cm
# -- days_to_flower = time from planting to flowering in days

# We can initially look at the effect of days to flower and treatment on plant size 

colors = ifelse(d$treatment == "low_water", "red", "blue") 
plot(stem_length~days_to_flower, d, pch=16, col=colors)

# Relationships look pretty non-linear. What if we log-transform the response variable? 
plot(log(stem_length)~days_to_flower, d, pch=16, col=colors)

m1 <- lm(log(stem_length)~days_to_flower*treatment, d)
summary(m1)

# Or you might prefer the display() function in Gelman & Hill's arm library which is much cleaner:
display(m1)

# Here's a quick way to check two assumptions of the normal linear model:
par(mfrow=c(2, 1), mar=rep(3,4), mgp=c(2,1,0))
plot(m1, which=1:2)
# Note: If you want to see an example of what these plots look like when these assumptions are badly violated, try refitting the model on non-log-transformed stem length and the repeating this plot. 

# But also note: per Gelman and Hill, these are really not terribly important assumptions and mild violations are no big deal! They don't even recommend looking at normal quantile plots of residuals. Residuals can however be useful to look at because they can reveal other things, such as weird outliers, nonlinearity, and grouping structure (i.e. non-independence) in the data. Much more on this later in the course.


#### Interpreting regression model parameters  ####

# How would you interpret the intercept and the effect of treatment in this model? # the intercept doesn't tell you very much except a ballpark estimate of where you'd expect the trend to be - obviously it can't take 0 days to flower
# the effect of treatment tells you that low water delays flowering, and limits stem length.

# Does this give an interpretation that makes practical sense? If not, what can we do? 
# This makes sense.

# One reason to center the explanatory variables on the mean is to make the effects more interpretable. 
d.center <- d
d.center$days_to_flower <- d$days_to_flower - mean(d$days_to_flower, na.rm=T)
# Or we can do the same think using a function from dplyr: 
# d.center <- mutate(d, days_to_flower = scale(days_to_flower, center=TRUE, scale=FALSE))

# Then we can repeat the same regression and compare. 
m2 <- lm(log(stem_length)~days_to_flower*treatment, d.center)
display(m1)
display(m2)
dev.off()
plot(log(stem_length)~days_to_flower, d.center, pch=16, col=colors)

# Which coefficients change? Now what do the intercept and treament effects mean? 
# the coefficients illustrate a deviation away from the mean value; a negative value means that the treatment effect (low water) decreases the days to flower relative to the mean

# Q how would you compare the sizes of the effect of treatment and days to flower? # treatment has a larger effect on stem length than days to flower

# G&H recommend centering all variables, and dividing continuous variables by 2x their standard deviation. This makes the effects of continuous variables more comparable to those of 1/0 categorical variables. 
# See G&H pp56-57.

# It is, however, much more common practice to just divide by one standard deviation (rather than 2), and is what people usually mean when they say the data were centered and scaled. It doesn't matter much, because either procedure puts the explanatory variables on the same scale, which is good for computation, and puts regression coefficients on a similar footing for comparison.

d.scaled <- d.center 
d.scaled$days_to_flower <- d.center$days_to_flower/(2*sd(d.center$days_to_flower))
# d.scaled <- mutate(d, days_to_flower = scale(days_to_flower, center=TRUE, scale=TRUE))

m3 <- lm(log(stem_length)~days_to_flower*treatment, d.scaled)
display(m3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Pause here for class discussion ####
# Please let the instructor know you have reached this point! 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#### PART 3: How can you tell if a model is good, or good enough?  ####

# a) Single-number summaries. 
# R^2 and other measures of fit or prediction 

# b) Fit 
# Plot fitted values vs observed.

# It can be helpful to plot the data with fitted lines. 

plot(log(stem_length)~days_to_flower, data=d.scaled, pch=16, col=colors)
betas <- coef(m3) # this extracts the coefficients from the linear model
abline(betas[1], betas[2], col="blue") # regression line for treatment=0
abline(betas[1]+betas[3], betas[2]+betas[4], col="red") # regression line for treatment=1


# Plotting the regression line plus confidence intervals.
# Keeping it simple, let's use just one explanatory variable.

m4 <- lm(log(stem_length)~days_to_flower, d.scaled)

# set the range of the explanatory variable to use for displaying the predictions
x.pred <- seq(min(d.scaled$days_to_flower)*1.1, max(d.scaled$days_to_flower)*1.1, by=0.05) # generates a sequence of values going just outside the range of observed x values

# The function predict() lets us predict the response variable from various levels of the explanatory variable. 
erodium.pred <- predict(m4, data.frame(days_to_flower=x.pred), se.fit=TRUE, interval="prediction") 
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

# Why are the confidence intervals for the regression line narrower than the confidence intervals for the model predictions? 
# Not really sure about this point

# Having done all that, a quicker way to plot confidence intervals is to use ggplot's statistical fit methods: 
ggplot(d.scaled, aes(x=days_to_flower, y=log(stem_length))) + geom_point() + stat_smooth(method=lm)


#### Simulating replicate data ####

# Can you simulate replicate data using the model? Remember what the regression model actually is, as a probability distribution, and what the parameters of that probability distribution are. 

# Hint 1: Note a short-cut is to use the function sim() in the arm library to generate random values of the fitted regression coefficients. You can then use this information to produce replicate simulated data sets. 

# Hint 2: the function rnorm(n, mean, sd) generates n random normally distributed data points with mean = "mean" and standard deviation "sd".  

# Comment: the function simulate() is what you'd probably ACTUALLY use to do this, but for now, try doing it "by hand" using sim() and rnorm()

