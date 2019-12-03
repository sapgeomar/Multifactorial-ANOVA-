# Sahed Ahmed Palash, Biological Oceanography, GEOMAR
# Master Thesis, Data Analysis
# Multi-factorial ANOVA for the biomass of euphausiacea as the function of D/N, onshore/offshore and oxygen bins

# set the working directory
#getwd()
setwd("/home/sahed/Desktop/office/stats/anova")

# import important packages
library(ggplot2)
library(MASS)
library(car)
library(nlme)
library(mgcv)
library(lattice)

#importing dataset
MyData <- read.csv(file="/home/sahed/Desktop/office/stats/anova/df_ euphausiacea.csv", header=TRUE, sep="\t")
MyData$bins <- cut(MyData$o2, breaks=c(0,10,20,50,150,250))
MyData$bins<-as.factor(MyData$bins)
MyData$on_off<-as.factor(MyData$on_off)
MyData$D_N<-as.factor(MyData$D_N)
attach(MyData)

# create a dotchart just to see the distribution pattern
#dotchart(MyData$biomass, groups = factor(MyData$bins), xlab="Biomass [mgCm⁻³]")
dotchart(MyData$biomass, groups = factor(MyData$on_off), xlab="Biomass [mgCm⁻³]")
dotchart(MyData$biomass, groups = factor(MyData$D_N), xlab="Biomass [mgCm⁻³]")

# create a table for treatment and replicates level
table(MyData$bins, MyData$on_off, MyData$D_N)                                                   # unequal sample size and design is not good 
# but its not a problem for ANOVA but it will consume the test power####
# create a histogram on biomass
hist(MyData$biomass)

# create a boxplot on treatment level
boxplot(MyData$biomass~MyData$bins)
boxplot(MyData$biomass~MyData$on_off)
boxplot(MyData$biomass~MyData$D_N)

# design the distribution on sample and overall mean
plot.design(MyData$biomass~MyData$bins)
plot.design(MyData$biomass~MyData$on_off)
plot.design(MyData$biomass~MyData$D_N)

# create a lattice plot for biomass as a function (bin/onshore and offshore) and (bins/day and night)
bwplot(MyData$biomass~MyData$bins|MyData$on_off)
bwplot(MyData$biomass~MyData$bins|MyData$D_N)

# extracting information about treatment levels
means<- tapply(MyData$biomass, interaction(MyData$bins, MyData$on_off, MyData$D_N), mean, na.rm=T)
round(means,digits = 0)

# model the data into ANOVA
model1<-aov(MyData$biomass~MyData$bins*MyData$on_off*MyData$D_N)

# check the homogeneity of variance
plot(resid(model1)~fitted(model1))
abline(h=0,lwd=2,lty=2,col="black")

# test for the homogeneity of variance 
fligner.test(MyData$biomass~MyData$bins)                                          # high p values (homogeneous varience)
fligner.test(MyData$biomass~MyData$on_off)                                        # low p values (variences are not homogeneous)
fligner.test(MyData$biomass~MyData$D_N)                                           # low p values (varience are not homogeneous)
fligner.test(MyData$biomass~interaction(MyData$bins,MyData$on_off, MyData$D_N))   # low p values (varience are not homogeneous)

# normality of errors
hist(resid(model1))

# test for normality
shapiro.test(resid(model1))                                                       # low p value (data are not normally distributed)

# chek the influential data points
plot(cooks.distance(model1), type="h")

# model output (ANOVA)
summary(model1)                                                                   # significant 
capture.output(anova(model1), file="anovaTable_euphausiacea.txt")
# post-hoc test (Tukey HSD for the pairwise to see the means which are significantly different than eac other)
TukeyHSD(model1)
capture.output(TukeyHSD(model1), file = "tukeyTable_euphausiacea.txt")
