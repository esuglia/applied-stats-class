# Script for Class 7

# Model comparison, overfitting, model averaging

# PLS 298 F2019
# Andrew Latimer

library(arm)
library(lme4)
library(AICcmodavg)
library(BMA)
library(dplyr)
library(sp)
library(maps)
library(BMA)


#### Part 1 -- Model complexity for a real data set ####

# Here's a slightly more realistic example using the char height data from Class 5. 
char = read.csv("char_with_fake.csv")

# This data set adds a fake, randomly generated predictor that is just noise (rnorm).

# Why are we doing this? Recall that the penalty term in AIC is designed to offset the improvement in fit you would typically get if you added a randomly generated explanatory variable to your model -- one that has no "true" relationship to the response variable. 


# We can train the model on the data from one fire and predict to the second. 
char.fit <- filter(char, Fire==1) # training data
char.holdout <- filter(char, Fire==4) # prediction data

# Fit the models
char.m1 <- lm(CharHt~Slope*Treated, data=char.fit)
char.m2 <- lm(CharHt~Slope*Treated*fake, data=char.fit)


# Compare raw sum of squared error and penalized fit terms

# mean squared error of the model fit
mean(resid(char.m1)^2); mean(resid(char.m2)^2)

# AIC
AICc(char.m1); AICc(char.m2)
AIC(char.m1, char.m2) # out of curiosity how much different is the comparison using AIC? 

# Does BIC give the same model ranking? 
BIC(char.m1, char.m2) 


# Does the model with the better fit to the training data also predict better in a new data set? 
#   Here we can assess the fit using mean squared errors.
mean((predict(char.m1, newdata=char.holdout) - char.holdout$CharHt)^2)
mean((predict(char.m2, newdata=char.holdout) - char.holdout$CharHt)^2)
# Which model does better in this cross-validation test? 

# Visual display of model fits and predictions
par(mfrow=c(2,2))
plot(predict(char.m1)~char.fit$CharHt, main="Model 1 Fit")
abline(0,1)
plot(predict(char.m2)~char.fit$CharHt, main="Model 2 Fit")
abline(0,1)
plot(predict(char.m1, newdata=char.holdout)~char.holdout$CharHt, main="Model 1 Predictive performance")
abline(0,1)
plot(predict(char.m2, newdata=char.holdout)~char.holdout$CharHt, main="Model 2 Predictive performance")
abline(0,1)

# Another you could try: Generate new "fake" data again, many times, and see how often the different model evaluation criteria favor a model that includes the fake data.


#### Part 2 -- Model averaging ####

# The R libary BMA can do model averaging across models with different combinations of predictors, using a penalized model fit criterion to determine the model weights.


#### Data setup ####
fish <- read.csv("fish_abundance_climate.csv")

# Data on fish abundance in streams across mostly Eastern North America from http://rsos.royalsocietypublishing.org/content/3/6/160093 

# Just to get a look at the data, we can map the sample sites
fish.sp <- SpatialPointsDataFrame(fish, coords=fish[,c("Longitude", "Latitude")])
par(mfrow=c(1,1)); plot.new()
map(database="state")
points(fish.sp, pch=16, col=fish.sp$Species, cex=0.7)

# For this exercise, let's focus on just one species of fish. 
species <- "Campostoma_oligolepis"
fish.co <- filter(fish, Species==species)

names(fish.co) # There are 10 environmental variables in this data set, plus latitude and longitude. 
# Having a lot of potential explanatory variables is a very common situtation, and raises the question of how to choose which to include in the model. There's simply no one best way to do this. But here are some key steps to take. 

# Evaluate collinearity of the data -- the correlation of different explanatory variables with each other. 
x <- fish.co[, 2:14]
round(cor(x), 2)
corrplot(cor(x), color=T)


# Which variables are highly correlated with each other? 

# In general, it's helpful to throw out very highly correlated predictors at the outset. Here, the temperature variables are highly correlated, and so are 3 of the precipitation variables. So let's trim out all but one of each of these groups. Ideally this is done using some biological knowledge about which variables are most likely to be important and/or interpretable. An alternative which we won't go into here is to use some kind of dimension reduction like principal components analysis. 
fish.co <- fish.co[,c("log_rarefied_abundance", "Latitude", "Longitude", "Max_temp", "Mean_diurnal_range", "Annual_precip", "Precip_seasonality")]

# This leaves us with 6 potential explanatory variables. Simplistically, we can start by comparing a "full model" versus a minimal model.
# check for missing data
sum(!complete.cases(fish.co))

m0 <- lm(log_rarefied_abundance~1, fish.co)
mfull <- lm(log_rarefied_abundance~., fish.co)
summary(mfull)
AIC(m0, mfull)
BIC(m0, mfull)

# To proceed, we could manually or automatically (see function stepAIC() ) step through various combinations of variables. Alternatively, with this number of predictors, we could actually test all combinations and find a "best." Alternatively, we could look across many/all combinations and rank them, and use some combination of these to make inferences. This latter approach is model averaging. Two libraries that exist to automate this are MuMIn ("Multi-Model Inference") and bma (Bayesian Model Averaging).

#### Model averaging ####

# To use bma, it's a bit of work to set up. We need to assemble the predictors and response manually first. A key motivation for multi-model inference is to improve prediction. So for this exercise we'll split the data into one subset for model fitting and another subset for model validation. 
n <- nrow(fish.co) # sample size
holdout <- seq(1, n, by=5) # hold out 20% of data -- every 5th row
fish.fit <- fish.co[-holdout,]
fish.holdout <- fish.co[holdout,]

# Setting up the data for model fitting 
y <- fish.fit$log_rarefied_abundance
x.avg <- fish.fit[,c("Longitude", "Latitude", "Max_temp", "Mean_diurnal_range", "Annual_precip", "Precip_seasonality")]

# Data for model validation 
x2.avg = fish.holdout[,c("Longitude", "Latitude", "Max_temp", "Mean_diurnal_range", "Annual_precip", "Precip_seasonality")]

# Fit the models
fish.avg = bicreg(x = x.avg, y = y)

# Look at the component models and how they're weighted
summary(fish.avg) 

# Some things to look at here are: as the models get more complex, what happens to the goodness of fit (R2), penalized fit (BIC) and relative probability or weights (post prob)


#### Comparing models ####

# Use the fitted model average to predict to the holdout locations
fish.avg.pred = predict(fish.avg, newdata=x2.avg)

# For comparison, let's refit the null and "full" models to just the fitting data (non including the holdout data)
fish_m0 <- lm(log_rarefied_abundance~1, fish.fit)
fish_mfull <- lm(log_rarefied_abundance~., fish.fit)

# And we can also add an obviously over-fitted model that includes all 2-way interactions: 
fish_m2way <- lm(log_rarefied_abundance~.^2, fish.fit)

# Calculate mean squared error (MSE) for each model, which measures within-sample fit
mean((predict(fish_m0) - fish.fit$log_rarefied_abundance)^2)
mean((predict(fish_mfull) - fish.fit$log_rarefied_abundance)^2)
mean((predict(fish_m2way) - fish.fit$log_rarefied_abundance)^2)
mean((predict(fish.avg, newdata=fish.fit)$mean - fish.fit$log_rarefied_abundance)^2)

# Calculate mean squared predictive error (MSPE) for each model, which measures out-of-sample predictive performance
mean((predict(fish_m0, newdata=fish.holdout) - fish.holdout$log_rarefied_abundance)^2)
mean((predict(fish_mfull, newdata=fish.holdout) - fish.holdout$log_rarefied_abundance)^2)
mean((predict(fish_m2way, newdata=fish.holdout) - fish.holdout$log_rarefied_abundance)^2)
mean((fish.avg.pred$mean - fish.holdout$log_rarefied_abundance)^2)

# How well does these 4 models do in within-sample fit, and in out-of-sample prediction?
# Can you see any downside to using model averaging? Is it harder to report or interpret the results from model averaging? 


#### Plotting model fit #### 

# We can also visualize the in-sample and out-of-sample performance of these 3 models
par(mfrow=c(2, 3), mar=c(2,2,2,2))
  # Display the fits of the three models 
plot(predict(fish_m0)~fish.fit$log_rarefied_abundance, main="Model 0 Fit", ylim=c(-0.1, 2))
abline(0,1)
plot(predict(fish_m2way)~fish.fit$log_rarefied_abundance, main="Complex Model Fit", ylim=c(-0.1, 2))
abline(0,1)
plot(predict(fish.avg, newdata=x.avg)$mean~fish.fit$log_rarefied_abundance, main="Model Avg Fit", ylim=c(-0.1, 2))
abline(0,1)
  # display the predictive or out-of-sample performance of the three models. 
plot(predict(fish_m0, newdata=fish.holdout)~fish.holdout$log_rarefied_abundance, main="Model 0 Prediction", ylim=c(-0.1, 2))
abline(0,1)
plot(predict(fish_mfull, newdata=fish.holdout)~fish.holdout$log_rarefied_abundance, main="Complex Model Prediction", ylim=c(-0.1, 2))
abline(0,1)
plot(predict(fish.avg, newdata=x2.avg)$mean~fish.holdout$log_rarefied_abundance, main= "Model Avg Prediction", ylim=c(-0.1, 2))
abline(0,1)
