## Elena Suglia
## Homework 1
## PLS 298, Applied stats modeling

library(tidyverse)
library(dplyr)
library(lattice)
library(vioplot)
library(lme4)
library(arm)
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

# Fit models
m1 = lm(uptake~logconc*Type, d)
summary(m1)
display(m1)


