plot_and_lmer <- function(x,dv_name) {
  colors<-wes_palette("Darjeeling1",43,type=("continuous"))
  colors[44]<-'#000000'

  p1 <- ggplot(x, aes(x = Trial, y = dv, colour = pair)) +
    geom_point(size = 2, alpha=.8) +
    geom_line(size = 1, alpha=.4) +
    geom_line(aes(y = fit), size = 1, alpha=.7) +
    geom_line(aes(y = fit_fixef), size = 2, alpha=.9, colour = "black") +
    stat_summary(fun = "mean", geom = "line", alpha=.5, colour = "black", size=2) +
    theme_classic() +
    scale_colour_manual(values = colors) +
    labs(y = dv_name) +
    labs(x = "Trial")
    #theme(legend.position = "top") +

  colors2<-wes_palette("Zissou1",2,type=("continuous"))
  p2 <- ggplot(x, aes(x=Trial, y= dv, colour=as.factor(Condition))) +
    stat_summary(fun="mean", geom="line", alpha=.5, size=1) +
    stat_summary(fun="mean", geom="point", alpha=.5, size=2) +
    #stat_summary(data=x, aes(x=Trial, y= dv, colour=as.factor(Condition)), geom="errorbar", fun.data=mean_se, alpha=1) + 
    stat_summary(data=x, aes(x=Trial, y= dv, colour=as.factor(Condition)), geom="ribbon", fun.data=mean_se, alpha=.2) + 
    # mean_cl_boot
    # geom_line(data=x, aes(x=Trial, y=fit_fixef, colour=as.factor(Condition)), linetype=2, size=1, alpha=.4) +
    geom_line(aes(x=Trial, y = fit_fixef, colour=as.factor(Condition)), linetype = 2, size = 2, alpha=.9) +
    theme_classic() +
    scale_colour_manual(values = colors2) +
    labs(colour = "Expressive") +
    theme(legend.position="top") + 
    # legend.title=element_text()
    labs(y = dv_name) +
    labs(x = "Trial")

  return(list(p1,p2))
}