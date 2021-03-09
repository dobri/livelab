#mada cc
library(reshape2)
library(ez)
library(TOSTER)
library(lsr)
library(ggplot2)
library(lme4)
library(Rmisc)
library(ggplot2)
library(psych)
library(wesanderson)
library(lmerTest)

##PREP
setwd('/Users/emilywood/desktop/r/madawaska')
imageDirectory<-"/Users/emilywood/desktop/r/madawaska'"
options(max.print=1000000)


##READ IN CC DATA
#position cc data
cc_pos<-read.csv("mada_cc_position.csv", header = TRUE,na.strings = c("", "NA"))
#detrended position cc data
cc_pos_dt<-read.csv("mada_cc_position_detrended.csv", header = TRUE,na.strings = c("", "NA"))
#acceleration cc data
cc_a<-read.csv("mada_cc_A.csv", header = TRUE,na.strings = c("", "NA"))

head(cc_pos)
head(cc_pos_dt)
head(cc_a)


#read in rating/averages data
avgs<-read.csv("mada_avgs.csv", header = TRUE,na.strings = c("", "NA"))
head(avgs)

#Correlations between average group gc and wcc and  performance ratings. 
#Check data
boxplot(avgs$wcc)
boxplot(avgs$gc)
boxplot(avgs$quality)
#Using Spearman because the performance quality ratings look wonky

cor.test(avgs$wcc,avgs$gc,method="spearman",exact=FALSE)
#sig negative relationship between average wcc and average gc

cor.test(avgs$wcc,avgs$quality,method="spearman",exact=FALSE)
#sig positive relationship between average wcc and performers' average ratings of performance quality
plot(avgs$wcc,avgs$quality)

cor.test(avgs$gc,avgs$quality,method="spearman",exact=FALSE)
#no relationship between average gc and performance quality ratings


##FACTORING
cc_pos$condition = factor(cc_pos$condition)
cc_pos$piece<-factor(cc_pos$piece, ordered=TRUE)
cc_pos$pair<-factor(cc_pos$pair, ordered=TRUE)

cc_pos_dt$condition = factor(cc_pos_dt$condition,labels=c("Mechanical","Expressive"))
cc_pos_dt$piece<-factor(cc_pos_dt$piece, ordered=TRUE)
cc_pos_dt$pair<-factor(cc_pos_dt$pair, ordered=TRUE)

cc_a$condition = factor(cc_a$condition)
cc_a$piece<-factor(cc_a$piece, ordered=TRUE)
cc_a$pair<-factor(cc_a$pair, ordered=TRUE)



#First, let's analyse acceleration, lines 70-180ish
#analysis of position starts around 180
#Dobri's plots showing individual slopes
g <- vector("list",4)
for (p in seq(1,2)) { #piece 1 and 2
  for (c in seq(1,2)) { #condition 1 and 2
    df<-cc_a[(cc_a$piece==p) & (cc_a$condition==c),]
    #colors<-wes_palette("FantasticFox1",length(unique(df$pair)),type=("continuous"))
    g[[p+2*(c-1)]] <- ggplot(df) + #plot gc against trial
      geom_jitter(aes(x=trial, y=wcc, colour=pair), size=2, alpha=.8, height=.00, width=.0) +
      geom_line(aes(x=trial, y=wcc, colour=pair), size=1.2, alpha=.5) +
      #scale_colour_manual(values=colors) +
      stat_summary(aes(x=trial,y=wcc), fun.y='mean', geom='line', size=1.2, alpha=.7) +
      stat_summary(aes(x=trial,y=wcc), geom="ribbon", fun.data=mean_cl_boot, alpha=.5) +
      #stat_summary(aes(x=trial,y=cc_a), geom="ribbon", fun.data=mean_se, alpha=.5) +
      theme_classic() +
      theme(axis.title = element_text(size=18),
            axis.text = element_text(size=14, colour="black")) +
      labs(x = "Trial", y = "Correlational-coupling", color="Musician pairing")+
      if (c == 1) {
        if (p==1){
          scale_x_continuous(breaks = c(1, 3, 5,7), labels=c("1","2","3","4"))
        } else {
          scale_x_continuous(breaks = c(2, 4, 6,8), labels=c("1","2","3","4"))
        }
      } else if (c == 2) {
        if (p==1){
          scale_x_continuous(breaks = c(2, 4, 6,8), labels=c("1","2","3","4"))
        }else{
          scale_x_continuous(breaks = c(1, 3, 5,7),labels=c("1","2","3","4"))
        }
      }
    head(df)
  }
}

print(g[1]) #p1,c1 (m)
print(g[2]) #p2,c1 (m)
print(g[3]) #p1,c2 (e)
print(g[4]) #p2,c2 (e)



#more plots, just showing average line

mada1<-cc_a[which(cc_a$piece==1),]
mada2<-cc_a[which(cc_a$piece==2),]

(cc_p1_a<-ggplot(data=mada1, aes(x=trial, y=wcc, color=condition)) +
    stat_summary(fun.y = mean, geom = "point", size=4)+
    stat_summary(fun.y = mean, geom = "line", aes(group = condition),size=1.25) + 
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + 
    scale_color_manual(values=c('firebrick4', 'goldenrod1'))+
    labs(x = "Trial", y = "Correlational-coupling (r)", color="Playing style")+
    #ggtitle("Effect of playing style over time on overall \ngranger-coupling (causal density)")+
    theme_bw()+
    scale_x_continuous(limits=c(1, 8), breaks=seq(1,8,1))+
    #scale_y_continuous(limits=c(0.005, 0.015)) +
    theme(panel.grid.major.x = element_line(),
          panel.grid.major.y = element_line(),
          plot.title = element_text(size = 22, vjust = 1.5),
          axis.title = element_text( size=18),
          legend.title = element_text(size = 18),
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size=16),
          legend.background = element_rect(color = "black"),
          axis.title.y = element_text(vjust= 1.8),
          axis.title.x = element_text(vjust= -0.5),
          axis.text = element_text(size=16, colour="black")))

#ggsave("cc_p1.png", width=14, height=11)

#same for piece 2
(cc_p2_a<-ggplot(data=mada2, aes(x=trial, y=wcc, color=condition)) +
    stat_summary(fun.y = mean, geom = "point", size=4)+
    stat_summary(fun.y = mean, geom = "line", aes(group = condition),size=1.25) + 
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + 
    scale_color_manual(values=c('firebrick4', 'goldenrod1'))+
    labs(x = "Trial", y = "Correlational-coupling (r)", color="Playing style")+
    #ggtitle("Effect of playing style over time on overall \ngranger-coupling (causal density)")+
    theme_bw()+
    scale_x_continuous(limits=c(1, 8), breaks=seq(1,8,1))+
    #scale_y_continuous(limits=c(0.005, 0.015)) +
    theme(panel.grid.major.x = element_line(),
          panel.grid.major.y = element_line(),
          plot.title = element_text(size = 22, vjust = 1.5),
          axis.title = element_text( size=18),
          legend.title = element_text(size = 18),
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size=16),
          legend.background = element_rect(color = "black"),
          axis.title.y = element_text(vjust= 1.8),
          axis.title.x = element_text(vjust= -0.5),
          axis.text = element_text(size=16, colour="black")))

#ggsave("cc_p2.png", width=14, height=11)

#LME for acceleration
m00=lmer(wcc ~ 1 + (1|pair),data=cc_a,REML=0)
m01=lmer(wcc ~ 1 + trial + (1|pair),data=cc_a,REML=0)
m02=lmer(wcc ~ 1 + trial + piece + (1|pair),data=cc_a,REML=0) 
m03=lmer(wcc ~ 1 + trial+condition + piece + (1|pair),data=cc_a,REML=0) 
m04=lmer(wcc ~ 1 + trial*condition + piece + (1|pair),data=cc_a,REML=0) 
m05=lmer(wcc ~ 1 + trial*condition*piece + (1|pair),data=cc_a,REML=0)  

anova(m00,m01,m02,m03,m04,m05)

#uh oh.. I get an error for each model 'boundary (singular) fit'
#I don't know how to deal with this so will come back to it later.



#Let's look at the position data
#Dobri's plots
g <- vector("list",4)
for (p in seq(1,2)) { #piece 1 and 2
  for (c in seq(1,2)) { #condition 1 and 2
    df<-cc_pos[(cc_pos$piece==p) & (cc_pos$condition==c),]
    #colors<-wes_palette("FantasticFox1",length(unique(df$pair)),type=("continuous"))
    g[[p+2*(c-1)]] <- ggplot(df) + #plot gc against trial
      geom_jitter(aes(x=trial, y=wcc, colour=pair), size=2, alpha=.8, height=.00, width=.0) +
      geom_line(aes(x=trial, y=wcc, colour=pair), size=1.2, alpha=.5) +
      #scale_colour_manual(values=colors) +
      stat_summary(aes(x=trial,y=wcc), fun.y='mean', geom='line', size=1.2, alpha=.7) +
      stat_summary(aes(x=trial,y=wcc), geom="ribbon", fun.data=mean_cl_boot, alpha=.5) +
      #stat_summary(aes(x=trial,y=cc_pos), geom="ribbon", fun.data=mean_se, alpha=.5) +
      theme_classic() +
      theme(axis.title = element_text(size=18),
            axis.text = element_text(size=14, colour="black")) +
      labs(x = "Trial", y = "Correlational-coupling", color="Musician pairing")+
      if (c == 1) {
        if (p==1){
          scale_x_continuous(breaks = c(1, 3, 5,7), labels=c("1","2","3","4"))
        } else {
          scale_x_continuous(breaks = c(2, 4, 6,8), labels=c("1","2","3","4"))
        }
      } else if (c == 2) {
        if (p==1){
          scale_x_continuous(breaks = c(2, 4, 6,8), labels=c("1","2","3","4"))
        }else{
          scale_x_continuous(breaks = c(1, 3, 5,7),labels=c("1","2","3","4"))
        }
      }
    head(df)
  }
}

print(g[1]) #p1,c1 (m)
print(g[2]) #p2,c1 (m)
print(g[3]) #p1,c2 (e)
print(g[4]) #p2,c2 (e)



#PLOT average lines

mada1<-cc_pos[which(cc_pos$piece==1),]
mada2<-cc_pos[which(cc_pos$piece==2),]

(cc_p1_pos<-ggplot(data=mada1, aes(x=trial, y=wcc, color=condition)) +
    stat_summary(fun.y = mean, geom = "point", size=4)+
    stat_summary(fun.y = mean, geom = "line", aes(group = condition),size=1.25) + 
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + 
    scale_color_manual(values=c('firebrick4', 'goldenrod1'))+
    labs(x = "Trial", y = "Correlational-coupling (r)", color="Playing style")+
    #ggtitle("Effect of playing style over time on overall \ncorrelational-coupling")+
    theme_bw()+
    scale_x_continuous(limits=c(1, 8), breaks=seq(1,8,1))+
    #scale_y_continuous(limits=c(0.005, 0.015)) +
    theme(panel.grid.major.x = element_line(),
          panel.grid.major.y = element_line(),
          plot.title = element_text(size = 22, vjust = 1.5),
          axis.title = element_text( size=18),
          legend.title = element_text(size = 18),
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size=16),
          legend.background = element_rect(color = "black"),
          axis.title.y = element_text(vjust= 1.8),
          axis.title.x = element_text(vjust= -0.5),
          axis.text = element_text(size=16, colour="black")))

#ggsave("cc_p1.png", width=14, height=11)


#same for piece 2
(cc_p2_pos<-ggplot(data=mada2, aes(x=trial, y=wcc, color=condition)) +
    stat_summary(fun.y = mean, geom = "point", size=4)+
    stat_summary(fun.y = mean, geom = "line", aes(group = condition),size=1.25) + 
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + 
    scale_color_manual(values=c('firebrick4', 'goldenrod1'))+
    labs(x = "Trial", y = "Correlational-coupling (r)", color="Playing style")+
    #ggtitle("Effect of playing style over time on overall \ngranger-coupling (causal density)")+
    theme_bw()+
    scale_x_continuous(limits=c(1, 8), breaks=seq(1,8,1))+
    #scale_y_continuous(limits=c(0.005, 0.015)) +
    theme(panel.grid.major.x = element_line(),
          panel.grid.major.y = element_line(),
          plot.title = element_text(size = 22, vjust = 1.5),
          axis.title = element_text( size=18),
          legend.title = element_text(size = 18),
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size=16),
          legend.background = element_rect(color = "black"),
          axis.title.y = element_text(vjust= 1.8),
          axis.title.x = element_text(vjust= -0.5),
          axis.text = element_text(size=16, colour="black")))

#ggsave("cc_p2.png", width=14, height=11)

#stats
m00=lmer(wcc ~ 1 + (1|pair),data=cc_pos,REML=0)
m01=lmer(wcc ~ 1 + trial + (1|pair),data=cc_pos,REML=0)
m02=lmer(wcc ~ 1 + trial + (1+trial|pair),data=cc_pos,REML=0)
#get boundary fit problem, so I'm going to drop  trial|pair
m02=lmer(wcc ~ 1 + trial + piece + (1|pair),data=cc_pos,REML=0) 
m03=lmer(wcc ~ 1 + trial+condition + piece + (1|pair),data=cc_pos,REML=0) 
m04=lmer(wcc ~ 1 + trial*condition + piece + (1|pair),data=cc_pos,REML=0) 
m05=lmer(wcc ~ 1 + trial*condition*piece + (1|pair),data=cc_pos,REML=0)  


anova(m00,m01,m02,m03,m04,m05)

summary(m02)
#wcc significantly increases across trial b=.002 
#wcc is significantly higher in piece 2 than piece 1

#check linearity
plot(resid(m02),cc_pos$wcc)
#doesn't look like it passes lienarity...
#Shouldn't these look more random? So can I even run this LME? I think I'm violating an assumption here. 

#check homogeneity
plot(m02)

#check residuals
require("lattice")
qqmath(m02)
qqnorm(resid(m02))



