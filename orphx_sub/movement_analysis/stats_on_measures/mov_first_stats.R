library(lme4)
library(texreg)
library(lmerTest)
library(ggplot2)
library(wesanderson)

# Point the working directory to where you've got the R files and the pre-processed data.
# setwd() 

# This does the stats and draws and saves figures.
source("plot_and_lmer.R")

X<-read.csv('rez_2021-06-20.csv')
X$Condition<-X$Condition-1
X$pp<-factor(X$pp)
summary(X)

sink(paste('diary_lmems_all_dvs_',Sys.Date(),'.txt',sep='')) # This will print all raw stats output to a text file.
for (dv in seq(5,13)){
  print(names(X)[dv])
  print(names(X)[dv])
  print(names(X)[dv])
  # Set the last arg to 1 to save raw figs with individual data. 2 for nicer summary plots with Trial and Condition.
  p<-plot_and_lmer(X,dv,0)
}
sink()
