#### PLS 298 Applied statistical modeling ####
#### Script file for Class session #1 -- Exploratory Data Analysis ####

# Note the script assumes that you have the class 1 data set saved into the same directory as the script.

# If you'd like to quickly toggle to different parts of this script, you can open the table of contents subwindow (button at top right of this window in RStudio)

install.packages("MASS")
install.packages("lattice")
install.packages("vioplot")

library(MASS)
library(lattice)
library(vioplot)
library(dplyr)
library(ggplot2)
library(dplyr)

# If you need to use an R library or package and don't already have it, you can tell R to download and install it. For example:
# install.packages("vioplot")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### PART 1 -- Some basic and useful tools for data visualization in R ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# To introduce some basic but very useful tools for visualizing data, we can start out using the Sitka data set from the MASS library. This data set is from an experiment in which Sitka spruce trees were grown with and without exposure to ozone. It contains the following:
# size -- height of the tree in meters
# Time -- the time the tree height was measured (there are >1 times per tree)
# tree -- the index number for each individual tree
# treat -- a factor variable indicating which treatment the tree received, control or ozone

# Check the data set is loaded and looks OK
head(Sitka)
class(Sitka) # tells you what kind of object this is. In this case, we want it to be a data frame (or a "tibble," which is a form of data frame)
str(Sitka) # more info about the object, including about the different columns it includes, and the types of date in each column. 

# First it's useful to look at the overall distribution of the data.
# That's good to know and also helps check for data entry errors.
hist(Sitka$size)
stem(Sitka$size)

# For example, what if someone typed "61.5" instead of "6.15"?
# Insert the error into the data set:
Sitka$size[which(Sitka$size==6.15)] <- 61.5
# Now it's easy to see it if we just look:
stem(Sitka$size)
# And if we're sure it's a typo we could fix it, or we could set the observation to NA (no data). 
Sitka$size[which(Sitka$size==61.5)] <- 6.15 # <- NA
stem(Sitka$size)

# Make quantile plots of the response variable. A normal quantile plot visually compares the distribution of some data to a normal distribution. 
qqnorm(Sitka$size, main="Normal Q-Q plot", xlab="Theoretical quantiles", ylab="Sample quantiles", )
	# What would it mean if points in the upper tail of the distribution were far above the trend line? 
  # A: skewed to the right
	# For example if the plot looked like this?
Sitka_size_skewed <- c(Sitka$size, rnorm(10, mean=8, sd=1))
qqnorm(Sitka_size_skewed, main="Normal Q-Q plot", xlab="Theoretical quantiles", ylab="Sample quantiles")
# Hint: hist(Sitka_size_skewed)
# Note quantile plots can be also be used to compare data to other distributions, or to compare two data sets.  

#### A) Looking at how the response variable Size is distributed within and across various subgroups of the data. ####

# We can manually break the data into groups for display. Usually we would use automated ways of doing this (see trellis plots and box plots below). But for now, just to make it very explicit, we'll do it manually.

# We are often interested in what the data look like in different experimental treatments. To do this, we can subset the data by treatment level and look at those subgroups. Here treatment is in the column "treat" and has two levels, which we can identify this way:
levels(Sitka$treat) 

x1 <- filter(Sitka, treat=="control")
  # Note the function filter::dplyr is a cleaner way to select rows of a data frame than the more direct but clunky equivalent of selecting elements directly: Sitka[Sitka$treat=="control",]
x2 <- filter(Sitka, treat=="ozone")
vioplot(x1$size, x2$size, names=c("control", "ozone"), col="lightgray")

# We can use the ability of ggplot to split the data into groups before plotting to make a similar plot: 
ggplot(Sitka, aes(y=size, x=treat)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), aes(fill=treat))

# Boxplots are a simple and widely used kind of visual summary and can be used to compare data across groups: 
boxplot(size~Time, Sitka, xlab="Time (days)")
	# Q: What are the different parts of the boxplot displaying and how might those be useful?
  # A: "the five-number summary": minimum, first quartile, median, third quartile, and maximum. The empty circles are outliers.
# Adding "notches" to a boxplot gives a rough visual way to see whether the means of the groups are significantly different: 
boxplot(size~Time, Sitka, xlab="Time (days)", notch=T)

# We can also make boxplots using the formula interface
boxplot(size~treat, Sitka, ylab ="Size (m)", xlab="Treatment", notch=T, col="lightgreen")
boxplot(size~treat+Time,Sitka, col=c("lightblue", "pink"), ylab="Size (m)")
# A quick note on display -- here the axis labels won't fit, so you can rotate them if you want
boxplot(size~treat+Time,Sitka, col=c("lightblue", "pink"), ylab="Size (m)", las=2)


#### B) Looking at relationships among variables in the data set ####

# Basic scatter plot
plot(size~Time, Sitka, xlab="Time (days)", ylab="Size (m)") 

# For looking at continuous data against a factor or grouping variable, in addition to boxplots there are stripplots (also called dotplots):
stripplot(Time~size, data=Sitka, jitter.data=T, group=treat, aspect=1, xlab="Size (m)")

# Lattice plots for displaying scatterplots of two factors stratified by a third
	# Q: What do these show you that the simple scatterplot doesn't? 
stripplot(Time~size|treat, data=Sitka, jitter.data=T, aspect=1, xlab="Size (m)")

# It's often useful to add a smoothed trendline for visualization (Chambers et al. recommend this -- do you think it helps?)
xyplot(size~Time|treat, Sitka, type=c("p", "smooth"))

# Note that a similar plot can also be made with the ggplot2 library.
plot1 <- ggplot(Sitka, aes(Time, size)) + geom_point() + facet_wrap(~treat) 
plot1
plot1 + geom_smooth(span=2)
# By changing the "theme", we can make a visually cleaner version of the same thing: 
plot1 + geom_smooth(span=2, se=FALSE) + theme_classic()

# Going further, there are a bunch of "geom" functions that add things to plots, like smoothed lines, statistical fits with/without error bars etc. 
# And while ggplot2 is a complex library, it allows changing all kinds of things about the appearance of the plot. 
plot1 + theme(panel.background=element_blank()) + theme(title=element_text(colour="purple", size=12))+ theme(axis.line = element_line(size=0.5, colour="orange3")) + theme(strip.background= element_rect(fill="pink")) + theme(strip.text=element_text(size=14, colour="green")) # etc. 
# Ready for publication! :) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### PART 2 -- In-class exercise:                      #### 
#			    Relationships among continuous variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here is another data set. At a few dozen grassland sites around the world, researchers measured aboveground annual net primary producticity (ANPP1). We also have the mean annual temperature and precipitation from these sites, as well as the type of ecoregion of each site. The data were gathered together to explore how productivity relates to climate.

# The data are from here: 
# http://dx.doi.org/10.3334/ORNLDAAC/654 

d <- read.csv("Grassland_NPP_31_site_summary.csv", na.strings = "-9999")
d <- select(d, Site, Name, Rainfall, Temperature, Ecoregion, ANPP1) # keep only the 

# Each row in the data set represents one site where ANPP was measured. At each site, there are some environmental observations (mean annual "Temperature" in degrees C, and "Rainfall" in mm). The broad Ecoregion is also listed. 

# We would like to use these data to learn something about the relationship between temperature and precipitation on the one hand, and primary productivity on the other. Use graphical techniques to explore to what extent we might be able to do this, and to identify weaknesses. 

# Some particular questions:

# 1) What is the distribution of ANPP across the sites? Does it appear to be normally distributed?

ggplot(d, aes(Site, ANPP1)) +
  geom_point()

hist(d$ANPP1)

qqnorm(d$ANPP1) 
# does not appear to be normally distributed
# there is one very large outlier ANPP1 at site krk that's skewing the data to the right; it may be better to remove it

# 2) Will we be able to make distict conclusions about whether precipition and temperature are associated with ANPP? if not, why not? What might improve our ability to do that? 

ggplot(d, aes(Temperature, ANPP1)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Ecoregion)

ggplot(d, aes(Rainfall, ANPP1)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Ecoregion)

# there appears to be a bit of a positive relationship between both temp & rainfall and ANPP1, but it is obscured a little by the outlier point

# 3) Does it appear that ANPP relates differently to rainfall and temperature in different types of ecosystems ("Ecoregion")? If so, what would be some ways to what extent can we take that variation into account, or nonetheless make inferences about ANPP responses to  precipitation and temperature?

# separating by Ecoregion breaks up the data into too few data points per Ecoregion to make robust conclusions
# however, in the drier regions it seems like rainfall is a limiting factor for NPP, and increased rain = increased NPP
# in the cold steppe if you squint hard enough there is higher NPP at higher temps; same with savanna



