# Script for Class 5
# PLS 298 F2019
# Andrew Latimer

library(arm) # library associated with Gelman & Hill book 
library(lme4)
library(lattice)
library(dplyr)

### Part 1: Grouping structure in data ####
#          Comparing complete pooling and no-pooling models

# This examples uses the radon data set from Gelman & Hill. This data set contains radon measurements from houses in 85 counties in Minnesota.
# The variables in the data set are: 
# y = log-transformed radon measurement
# x = where the measurement was taken (0 = basement, 1 = ground floor)
# county = county in which house is located

radon <- read.csv("radon.csv")

# A couple things to take a look at the data
hist(radon$y)
boxplot(y~county, radon)
stripplot(county~y, data=radon, scales=list(cex=0.5))
# what does this tell us about where the variation is (within counties? among counties?). What does it reveal about the sample sizes?  
# more variation within than between counties
# there are some counties that have far greater sample sizes

# To look at the distribution of sample sizes for all the counties, we could do this:
hist(tabulate(radon$county), col="ivory1")
# Or this as a more community-ecology oriented display:
plot(rev(sort(tabulate(radon$county))), ylab="Sample size", xlab="County ranked by sample size", pch=16, col="darkgoldenrod2")
# Is the sampling balanced? 
# no

# Let's compare alternative ways to analyze these data:
# a) complete pooling, no grouping variable (just ignore the fact there are groups)

m1 <- lm(y~x, radon)
display(m1)
# more radon in the basement (neg correl b/w radon meas & floor 0/1)

# We can check out the residuals to see if there's any obvious structure in them:
plot(resid(m1), pch=16, col="cadetblue4")
boxplot(resid(m1)~county, data=radon)
# It's hard to see because there are so many counties, but looks like there is county-related variation in the residuals. So we should take this into account somehow. 

# b) No-pooling model -- with a categorical variable for county.
#     This is not pooling because effects are estimated separately for each county. 

m2 <- lm(y~x+county, radon)
summary(m2)
display(m2)
# This is clearly a better-fitting model. Some of the county basement radon estimates (intercept plus county effect) are large and significant, though most are not.
plot(coef(m2))
# One thing to notice is that some of the biggest estimated county effects are from counties with small sample sizes. County "LAC QUI PARLE" with the largest coefficient has only 2 data points! 
grep("LAC QUI PARLE", radon$county)
# or maybe more usefully: 
sort(table(radon$county))
# So the county effects, though they might come out as statistically significant, aren't very confidence inspiring. Gelman and Hill comment that they are effectively overfitted. 

# One illustration of this is to compare effect sizes to sample sizes. Would you expect estimated effect sizes to be related to sample size? How and why? 
samplesize <- tabulate(radon$county)
m2.coef <- coef(m2)
m2.coef
effectsize <- c(0, m2.coef[3:86]) + m2.coef[1] # extract the effect sizes for the levels of categorical variable "county", tacking on a 0 because (in lm and lmer) the first county is by default used as the base level against which the other levels are contrasted. 
plot(samplesize, effectsize)
# This kind of plot is sometimes called a "funnel plot" because of the shape. In general, unusually large effect sizes tend to be associated with small sample sizes. This has caused problems for many fields in which sample sizes tend to be small, including notoriously psychology, but also ecology.


# c) Random effects model as a compromise between pooling and no-pooling estimates. 
m3 <- lmer(y~x+(1|county), data=radon)

# compare the effects for county in this model
ranefsize <- coef(m3)$county$`(Intercept)` # extract the coefficient values from the fitted mixed-effects model object
plot(effectsize, ylim=c(0, 3), pch=16, col="darkgray") # plot the coefficients for county from the fixed-effects model (the no-pooling model above) 
points(ranefsize, pch=16, col="orange")
abline(fixef(m3)[1], 0, lty=2) # add a line at the overall mean radon level
# What is the pattern here, when you compare the gray to the orange points? 

# Does the degree of "shrinkage" depend on sample size? 
effect.diff <- abs(effectsize-ranefsize)
plot(samplesize, effect.diff, ylab="amount of shrinkage", xlab="sample size", pch=16, col="slateblue")
# What does this plot show?




##### Part 2: Using random effects to account for grouping/nonindependence in data ####

char = read.table("charheightdatav2.txt", header=T)

# In this data set there are measurements of char heights measured on tree trunks after fire. The data are all from the same fire, but were collected at plots grouped into transects.In other words, the individual plots weren't placed at random around the fire. The locations of the transects was randomized.
# The variables are: 
# CharHt = Char height in meters, an indicator of fire intensity
# Transect = categorical transect ID
# Steepness = steepness of the slope of the terrain at each plot, centered and scaled
# Treated = whether the forest had been thinned before it burned (we can ignore this for now)

# For this example, let's look just at a single fire: 
char <- filter(char, Fire == 1) 

# Is there a relationship betwen char height and slope steepness? 
plot(CharHt~Steepness, char)
xyplot(CharHt~Steepness|Transect, char)

# If you fit a model without including information about Transects, do you see obvious structure in the residuals? 

# What random effects structure would you use to model the relationship between steepness and char height?

# Does adding random effects to account for grouping in the data change the estimate of the effect of steepness? How many standard errors the steepness parameter is away from 0? 

# How much of the variation in individual char height measurements 
#   corresponds to differences among transects, versus within transects?
#   See G&H page 258.

