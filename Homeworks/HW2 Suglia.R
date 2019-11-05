
library(arm)
library(lme4)
library(lattice)
library(tidyverse)
library(dplyr)
library(vioplot)
library(lme4)
library(arm)
library(bbmle) # Ben Bolker's library of mle functions
library(MASS)


# Question 1. 
# Using the Char Height data set from class (char_with_fake.csv), construct a model with random effects for these data, using CharHt as the response variable. Assume that Transect is the only relevant grouping factor, and that Steepness (of the topography) and Diameter (of the trees) are the only available predictors. In reporting about the model please include:
# •	a) A brief explanation of how you chose variables, and which (if any) you decided to allow to vary by group (Transect). 
# •	b) An assessment of how much variation there is in the group-level random effects. 
# •	c) A brief assessment of how well your selected model fits the data. 

# Notes about this data set:
  # This data set adds a fake, randomly generated predictor that is just noise (rnorm).

  # Why are we doing this? Recall that the penalty term in AIC is designed to offset the   improvement in fit you would typically get if you added a randomly generated explanator   y variable to your model -- one that has no "true" relationship to the response          variable. 

d = read.csv("char_with_fake.csv", header = TRUE)
head(d)

# Look at variation in data across groups
boxplot(CharHt~Transect, d, xlab = "Transect")
# Tons of variation among groups

# there appears to be some variation in response to both explanatory variables
ggplot(d, aes(x=Diameter, y=CharHt)) +
  geom_point() +
  geom_smooth(span=2) +
  theme_classic() +
  ggtitle("Diameter vs CharHt") +
  xlab("Diameter") +
  ylab("CharHt") +
  facet_wrap(~Transect)

ggplot(d, aes(x=Slope, y=CharHt)) +
  geom_point() +
  geom_smooth(span=2) +
  theme_classic() +
  ggtitle("Slope vs CharHt") +
  xlab("Slope") +
  ylab("CharHt") +
  facet_wrap(~Transect)

# Normality?
qqnorm(d$CharHt, main = "Normal Q-Q plot", xlab = "Theoretical quantiles", ylab = "Sample quantiles")
# seems ok

# What about collinearity?
x = select(d, Slope, Diameter)
round(cor(x), 2)
# They're only correlated a little (0.2); I think this is ok

# Let's fit some models and compare them

# First, decide which predictor(s) to include

# Fixed effects only
m1 = lm(CharHt~Slope, data =d)
m2 = lm(CharHt~Diameter, data = d)
m3 = lm(CharHt~Slope+Diameter, data = d)
m4 = lm(CharHt~Slope*Diameter, data = d)

# Mixed models including both fixed and random effects
m5 <- lmer(CharHt~Slope+(1|Transect), data=d, REML = FALSE)
m6 <- lmer(CharHt~Diameter+(1|Transect), data=d, REML = FALSE)
m7 <- lmer(CharHt~Diameter+(1+Slope|Transect), data=d, REML = FALSE)
m8 <- lmer(CharHt~Slope+(1+Diameter|Transect), data=d, REML = FALSE)
m9 <- lmer(CharHt~Diameter+Slope+(1|Transect), data=d, REML=FALSE)
m10 <- lmer(CharHt~Diameter*Slope+(1|Transect), data=d, REML = FALSE)

AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
BIC(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
# those that include the interaction are all worse than those without the interaction

# Stepwise check through all combinations of variables
# check for missing data
#sum(!complete.cases(d)) # no missing data

#m0 <- lm(CharHt~1, d)
#mfull <- lm(CharHt~., d)
#summary(mfull)
#AIC(m0, mfull)
#BIC(m0, mfull)

#stepAIC(d, )

# How much variation is there in the group-level random effects?
display(m8)

# Do we include the random effects or not?
# Notes for me from the help page:
# The marginal R2 value represents the variance explained by the fixed effects, defined as: 
# R_GLMM(m)² = (σ_f²) / (σ_f² + σ_α² + σ_ε²)
# The conditional R2 value is interpreted as the variance explained by the entire model, including both fixed and random effects, and is calculated according to the equation: 
# R_GLMM(c)² = (σ_f² + σ_α²) / (σ_f² + σ_α² + σ_ε²)
# where σ_f² is the variance of the fixed effect components, σ_α² is the variance of the random effects, and σ_ε² is the “observation-level” variance.

library(MuMIn)
r.squaredGLMM(m5) 
r.squaredGLMM(m6) 
r.squaredGLMM(m7) 
r.squaredGLMM(m8) 
r.squaredGLMM(m9) 
r.squaredGLMM(m10) 

# could use dredge to combine all combinations and rank them all

# Including the random effects does allow the model to explain more of the variation in the response variable.

# Testing model assumptions
par(mfrow=c(2, 1), mar=rep(3,4), mgp=c(2,1,0))
plot(m4, which=1:2)

# What about centering and scaling, in order to interpret the coefficients of the model? Should we write down our interpretations of the coefficients?

# Another way to assess model performance is with cross-validation
# The following example is outlined in Class 7 Script

# There are 6 fires in total: let's withhold 3 and allow the model to predict the last 3
d.fit <- filter(d, Fire %in% (1:3)) # training data
d.holdout <- filter(d, Fire %in% (4:6)) # prediction data

# Or, as in class script 7, use data on 1st fire to predict 2nd fire
d.fit <- filter(d, Fire ==1) # training data
d.holdout <- filter(d, Fire ==2) # prediction data

# Fit the models
char.m1 <- lmer(CharHt~Diameter+(1+Diameter|Transect), data=d.fit)
char.m2 <- lmer(CharHt~Diameter*fake+(1+Diameter|Transect), data=d.fit)

# Compare raw sum of squared error and penalized fit terms

# mean squared error of the model fit
mean(resid(char.m1)^2)
mean(resid(char.m2)^2)

mean((predict(char.m1, newdata=d.holdout, allow.new.levels = TRUE) - d.holdout$CharHt)^2)
mean((predict(char.m2, newdata=d.holdout, allow.new.levels = TRUE) - d.holdout$CharHt)^2)
# Which model does better in this cross-validation test? 



# Visual display of model fits and predictions
par(mfrow=c(2,2))
plot(predict(char.m1)~d.fit$CharHt, main="Model 1 Fit")
abline(0,1)
plot(predict(char.m2)~d.fit$CharHt, main="Model 2 Fit")
abline(0,1)
plot(predict(char.m1, newdata=d.holdout, allow.new.levels = TRUE)~d.holdout$CharHt, main="Model 1 Predictive performance")
abline(0,1)
plot(predict(char.m2, newdata=d.holdout, allow.new.levels = TRUE)~d.holdout$CharHt, main="Model 2 Predictive performance")
abline(0,1)

# I'm not sure what is happening here; it appears that both models are good. Perhaps the signal in the data is strong enough that the model can handle the noise produced by the randomly generated data?

# Another you could try: Generate new "fake" data again, many times, and see how often the different model evaluation criteria favor a model that includes the fake data.
# Maybe it's worth trying this...

# There is another example in Class 8 Script that withholds different parts of the data - the equivalent in this data set would be withholding half the transects...

# To set up this comparison, let's first break the radon data set into two parts by county, so we can make predictions to counties that are "new", that is, not contained in the fitting data set. 

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

# Model averaging?

## Question 2
# Get data on performance scores for pairs figure skating in the 1932 olympics (from http://www.stat.columbia.edu/~gelman/arm/examples/olympics/olympics1932.txt). This is formatted for R as “olympics.csv” on Smartsite in the homework folder.

# Let’s assume the question is: which is the bigger source of variation in the scores for skating programs, the judges or the skating pair? Fit a mixed model for these data with “score” as the response variable, and random effects for judge and skating pair (“judge” and “pair”). Interpret the results (coefficients and their standard errors, standard errors of the random effects). Is there a judge that tends to give consistently higher scores?

d2 = read.csv("olympics.csv", header = TRUE)
ggplot(data=d2, aes(x=judge, y=score, group = pair, color=pair)) +
  geom_point()
boxplot(score~judge, d2, xlab = "Judge")
# just by looking at the graphs, judge 7 appears to give consistently higher scores than the others, and perhaps judge 2 as well
boxplot(score~pair, d2, xlab = "pair")
# out of curiosity, looked at variation among pairs; not much variation of scores within pairs
# Fit model
m16 <- lmer(score~ (1|judge), data = d2)
mskate <- lmer(score~  (1|judge) + (1|pair), data = d2)
# Interpet results
summary(mskate)

#             coef.est coef.se
# (Intercept)  6.03     0.31  
# pair        -0.20     0.04  
# judge       -0.04     0.06  

# Error terms:
#   Groups   Name        Std.Dev.
# judge    (Intercept) 0.30    
# pair     (Intercept) 0.18    
# Residual             0.27

# standard error of coefficients for pair and judge are similar but slightly higher for judge
# pairs are ranked in descending order of average score; so it makes sense that there is a negative correlation between pair and score
# the number assigned to judge doesn't really mean anything so the coefficient for slope between judge and score is not easily interpretable

# the standard deviation in the error term for the random effect judge is again higher than that for pair
# I would conclude that judges explain more of the variation in scores than pairs do

