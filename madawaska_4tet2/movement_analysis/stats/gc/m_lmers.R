m_lmers <- function(x) {
  m<-list("vector",6)
  m[[1]]<-lmer(dv ~ 1+(1|pair), data=x, REML=0)
  m[[2]]<-lmer(dv ~ 1+Condition+(1|pair), data=x, REML=0)
  m[[3]]<-lmer(dv ~ 1+Condition+Trial+(1|pair), data=x, REML=0)
  m[[4]]<-lmer(dv ~ 1+Condition*Trial+(1|pair), data=x, REML=0)
  m[[5]]<-lmer(dv ~ 1+Condition*Trial+piece+(1|pair), data=x, REML=0)
  
  print(anova(m[[1]],m[[2]],m[[3]],m[[4]],m[[5]]))
  print(summary(m[[3]]))
  print(summary(m[[4]]))
  print(summary(m[[5]]))
  print(screenreg(list(m[[1]],m[[2]],m[[3]],m[[4]],m[[5]])))
  
  best_model=1
  f<-which((anova(m[[1]],m[[2]],m[[3]],m[[4]],m[[5]])[8]<.05))
  if (length(f)>0){
    best_model<-f[length(f)]  
  }
  print(sprintf('The top model is #%i.',best_model))
  
  fitted <- predict(m[[best_model]])
  fitted_fixef <- getME(m[[best_model]],'X') %*% fixef(m[[best_model]])
  
  return(list(fitted,fitted_fixef,m,best_model))
}