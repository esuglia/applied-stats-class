# Script Class 8
# Intro to hierarchical models part 1
# PLS 298 -- F2019
# Andrew Latimer

library(lme4)
library(arm)
library(ggplot2)
library(gridExtra)
library(dplyr)

#### Part 1: Fitting a simple 2-level hierarchical model with lmer() ####

radon <- read.csv("radon_cty.csv")
head(radon)

# This is the same radon data set as before, except now there is an extra column for county-level uranium content of the bedrock. 

# Because there were some counties with no uranium data, these were excluded, and now we have data for only 76 counties.

counties <- unique(radon$county) # list of counties in data set
n.counties <- length(counties)

# As before, we can make a random effects model in which variation among individual houses can be explained using the explanatory variable x and a random intercept for county. 

m1 <- lmer(y ~ x + (1|county), data=radon, REML=FALSE)
display(m1)

# Now we can also try to explain the variation among counties, using the county-level uranium measurements. (See Gelman & Hill page 266).

head(radon, 10) # the explanatory variable "u.full" has a value for every row (every individual house), but the value is the same for all houses within a given county. So it's really a group-level explanatory variable. 

# The model with the group-level variable added: 
m2 <- lmer(y ~ x + u.full + (1|county), data=radon, REML=FALSE)
display(m2)

#### Questions Part 1 ####
# Compare the county-level standard deviation term in the two models.

# a) Does including county-level bedrock uranium content reduce the county-level error term? Does it reduce the residual error term? Why or why not? 

# b) Does including county-level uranium levels improve the model? 



#### Part 2: Predicting from a multilevel model ####

# Predicting and simulating from a multilevel model are more complicated than for a single-level model, because there are multiple error terms in the model. This means you have to think about variation at more than one level when predicting or simulating. 

# For example, we might want to predict levels of radon at the county level, given some information about bedrock uranium content, or we might want to predict radon levels in individual houses within the county. If we have information at both individual house level and county level, then we will make these predictions differently.

# Let's say we want to predict radon levels in individual houses. How we predict these levels depends on whether the houses are in a county we already know something about, or in a new unknown county. 

# 1) If the houses are in counties already included in our data set for fitting the model, then we have random effect values for those counties. To predict house-level radon levels, we would use all the information from our model, including both fixed and random effects. We'd expect the predictions to be relatively good.

# To set up this comparison, let's first break the radon data set into two parts by county, so we can make predictions to counties that are "new", that is, not containes in the fitting data set. 

radon1 <- filter(radon, county %in% counties[1:50])
radon2 <- filter(radon, county %in% counties[51:n.counties])

# Then we can hold out 10% of data from the part of the data set that we will use to fit the model, so we can also predict to "known" counties -- counties that are in the part of the data set used to fit the model. 
holdout1 <- seq(2, as.integer(nrow(radon1)), by=10)
radon1.fit <- radon1[-holdout1,]
radon1.holdout <- radon1[holdout1,]
# Drop any holdout data that is from a county not included in the fitting data. 
radon1.holdout <- filter(radon1.holdout, county %in% unique(radon1.fit$county))

# Fit our model to the data from the "known" counties, minus the holdout data. This is the same as model m2, above, but fitted to a subset of the data
m.sub <- lmer(y ~ x + u.full + (1|county), data=radon1.fit, REML=FALSE)

fits.re <- predict(m.sub, newdata=radon1.holdout) # by default predict includes the random effects in predictions
pdata <- data.frame(y=radon1.holdout$y, fits.re)
p.full <- ggplot(pdata, aes(fits.re, y)) + geom_point() + geom_smooth(method=lm) + theme_bw()
p.full

# 2) If the houses are in new counties that are not included in the data we used to fit the model, then we have no random effect estimates. To predict radon levels in these houses, we could use only the covariate values and coefficients. 

fits.fixed <- predict(m.sub, newdata = radon2, re.form=NA) # tell predict() not to include random effects in prediction
pdata2 = data.frame(y = radon2$y, fits.fixed)
p.fixed <- ggplot(pdata2, aes(fits.fixed, y)) + geom_point() + geom_smooth(method=lm) + theme_bw()
grid.arrange(p.full, p.fixed, nrow=1, ncol=2)

#### Questions Part 2 ####
# Do the predictions for the known counties look better than the predictions for the new counties? How can you tell? 

#### Quantifying model fit in simple multilevel models ####

# While there is no simple R2 measure for these models, there is a rough equivalent. 
# (Nakagawa & Schlielzeth 2013, Methods in Ecol & Evol)
# As with the predictions above, to assess model fit, we have to decide whether to include the random effects or not. 
# If the main point of the R2 is to tell how well the explanatory variables contribute to explaining the observed variation, then it's important not to include the random effects. Nakagawa and Schielzeth call this "marginal R2". To get a sense of the overall model fit, then random effects can be included as well as well as fixed effects; this is called "conditional R2". 

library(MuMIn)
r.squaredGLMM(m2)



####Part 3 Crossed vs nested random effects examples ####

# These examples are adapted directly from Chapter 2 of Bates' book  "lme4: Mixed-effects Modeling with R"
# http://lme4.r-forge.r-project.org/book/Ch2.pdf

#### Crossed random effects ####

data(Penicillin)
head(Penicillin)

# These data contain the diameter of bacterial colonies (i.e. growth rate) from 6 bacterial samples (sample A-F), each grown on 24 replicate plates (plate a-x).

hist(Penicillin$diameter)

# In these data, each sample and each plate is used multiple times, but there is only one data point for each plate x sample combination. Therefore we can't fit nested or interacting random effects.  

xtabs(~ sample + plate, Penicillin)

# Instead, we can use crossed random effects to account separately for the effect of Plate and the effect of sample.

fm1 <- lmer(diameter ~ 1 + (1|sample) + (1|plate), data=Penicillin)
display(fm1)

# Is there more variation among plates or among samples? How much of the variance is explained by the model? 

#### Nested random effects ####

data(Pastes)
head(Pastes)
# This data set on the strength of an industrially produced chemical paste, made in large batches and then packed into individual casks. Two casks from each batch are randomly selected and tested for quality control. 
hist(Pastes$strength)
xtabs(~ batch + cask, Pastes)

# In this data set, each cask is unique to a batch -- is nested within it, in other words. There were 2 samples taken from each cask. Here, it doesn't make sense to use crossed random effects, because "cask a" means something different within each batch. So we need to fit nested random effects, which is essentially fitting random effects for the interaction of batch and cask. 

pastes.m1 <- lmer(strength ~ 1 + (1|batch/cask), Pastes, REML=F)
display(pastes.m1)



