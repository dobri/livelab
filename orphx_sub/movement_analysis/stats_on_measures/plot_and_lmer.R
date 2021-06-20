plot_and_lmer <- function(x,dv,saveflag) {
  colors<-wes_palette("Darjeeling1",43,type=("continuous"))
  colors[44]<-'#000000'

  x<-x[(is.na(x[,dv]))==0,]
  x$out<-x[,dv]
  iv <- names(x[dv])
  
  m1<-lmer(out ~ 1+(1|pp), data=x, REML=0)
  m2<-lmer(out ~ 1+Condition+(1|pp), data=x, REML=0)
  m3<-lmer(out ~ 1+Condition+Trial+(1|pp), data=x, REML=0)
  m4<-lmer(out ~ 1+Condition*Trial+(1|pp), data=x, REML=0)
  
  print(anova(m1,m2,m3,m4))
  print(summary(m3))
  print(summary(m4))
  print(screenreg(list(m1,m2,m3,m4)))
  
  x$fit <- predict(m3)
  x$fit_fixef <- getME(m3,'X') %*% fixef(m3)
    #fixef(m3)[1]+fixef(m3)[2]*x$Condition+fixef(m3)[3]*x$Trial
  p1 <- ggplot(x, aes(x = Trial, y = out, colour = pp)) +
    geom_point(size = 2, alpha=.8) +
    geom_line(aes(y = fit), size = 1, alpha=.7) +
    geom_line(aes(y = fit_fixef), size = 2, alpha=.9, colour = "black") +
    stat_summary(fun = "mean", geom="line", alpha=.5, colour="black", size=2) +
    theme_classic() +
    scale_colour_manual(values = colors) +
    labs(y = iv) +
    labs(x = "Trial")
    #theme(legend.position = "top") +
  if (saveflag==1){
    filename=paste("regs_ind_fit_and_averages_dv",dv,"_",Sys.Date(),'.png',sep='')
    ggsave(filename,width=4,height=3,dpi=300)
  }
  
  p2 <- ggplot(x, aes(x = Trial, y = out, colour = pp)) +
    geom_point(size = 2, alpha=.6) +
    geom_line(size = 1, alpha=.5) +
    stat_summary(fun = "mean", geom="line", alpha=.5, colour="black", size=2) +
    theme_classic() +
    scale_colour_manual(values = colors) +
    labs(y = iv) +
    labs(x = "Trial")
  if (saveflag==1){
    filename=paste("regs_ind_and_averages_dv",dv,"_",Sys.Date(),'.png',sep='')
    ggsave(filename,width=4,height=3,dpi=300)
  }
  
  colors2<-wes_palette("Zissou1",2,type=("continuous"))
  p3 <- ggplot(x, aes(x=Trial, y=out, colour=as.factor(Condition))) +
    stat_summary(fun="mean", geom="line", alpha=.5, size=1) +
    stat_summary(fun="mean", geom="point", alpha=.5, size=2) +
    #stat_summary(data=x, aes(x=Trial, y=out, colour=as.factor(Condition)), geom="errorbar", fun.data=mean_se, alpha=1) + 
    stat_summary(data=x, aes(x=Trial, y=out, colour=as.factor(Condition)), geom="ribbon", fun.data=mean_se, alpha=.2) + 
    # mean_cl_boot
    # geom_line(data=x, aes(x=Trial, y=fit_fixef, colour=as.factor(Condition)), linetype=2, size=1, alpha=.4) +
    theme_classic() +
    scale_colour_manual(values = colors2) +
    labs(colour = "VLF") +
    theme(legend.position="top") + 
    # legend.title=element_text()
    labs(y = iv) +
    labs(x = "Trial")
  if (saveflag==2){
    filename=paste("regs_averages_and_model_dv_",dv,"_",Sys.Date(),'.png',sep='')
    ggsave(filename,width=4,height=3,dpi=300)
  }
  
  return(list(p1,p2,p3))
}