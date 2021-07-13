library(lme4)
library(texreg)
library(lmerTest)
library(ggplot2)
library(wesanderson)

# Point the working directory to where you've got the R files and the pre-processed data.
setwd('~/logos/mcmc/livelab/madawaska_4tet2/movement_analysis/')

# This does the stats and draws and saves figures.
source('stats/gc/m_lmers.R')
source('stats/gc/plot_and_lmer.R')
source('stats/gc/multiplot.R')

# Import the data
x<-read.csv('data/excel/X_clean_processed_gc.csv')
xp<-read.csv('data/excel/Xpcs_processed_gc.csv')
xd<-read.csv('data/excel/X_detrended_processed_gc.csv')
v<-read.csv('data/excel/V_processed_gc.csv')

X<-x
X$gc_x<-x$gc
X$gc_pca<-xp$gc
X$gc_xd<-xd$gc
X$gc_v<-v$gc

cor(as.matrix(cbind(X$gc_x,X$gc_pca,X$gc_xd,X$gc_v)))

# """
# [,1]       [,2]       [,3]       [,4]
# [1,] 1.0000000 0.21232039 0.82863274 0.19309643
# [2,] 0.2123204 1.00000000 0.08425262 0.15739867
# [3,] 0.8286327 0.08425262 1.00000000 0.07005823
# [4,] 0.1930964 0.15739867 0.07005823 1.00000000
# """

X$Condition<-X$condition-1
X$Trial<-X$trial-1
X$pair<-factor(X$pair)
summary(X)

#sink(paste('diary_lmems_all_dvs_',Sys.Date(),'.txt',sep='')) # This will print all raw stats output to a text file.
for (dv in seq(9,12)){
  print(names(X)[dv])
  X$dv<-X[,dv]
  
  # Find the best model. Models in M[[3]]. The best one is indicated in M[[4]].
  M <- m_lmers(X)
  X$fit<-M[[1]]
  X$fit_fixef<-M[[2]]

  g_ind<-list("vector",2)
  g_ave<-list("vector",2)
  for (p in unique(X$piece)){
    print(paste('Piece',p))
    pp<-plot_and_lmer(X[X$piece==p,],names(X)[dv])
    g_ind[[p]]<-pp[1]
    g_ave[[p]]<-pp[2]
  }
  multiplot(plotlist=g_ind,layout=matrix(c(1,2),nrow=1,byrow=TRUE))
  multiplot(plotlist=g_ave,layout=matrix(c(1,2),nrow=1,byrow=TRUE))
  
  if (TRUE) {
    filename=paste("regs_ind_and_model_dv_",names(X)[dv],"_",Sys.Date(),'.png',sep='')
    png(filename=filename,width=8,height=6,units="in",res=300)
    multiplot(plotlist=g_ind,layout=matrix(c(1,2),nrow=1,byrow=TRUE))
    dev.off()
    
    filename=paste("regs_means_dv_",names(X)[dv],"_",Sys.Date(),'.png',sep='')
    png(filename=filename,width=8,height=6,units="in",res=300)
    multiplot(plotlist=g_ave,layout=matrix(c(1,2),nrow=1,byrow=TRUE))
    dev.off()
  }
}
