# setwd('~/logos/mcmc/livelab/orphx_sub/movement_analysis/stats_on_measures')
x<-c(0.5,0.458333333,0.416666667,0.527777778,0.486111111,0.5,0.583333333,0.569444444,0.5,0.486111111,0.444444444,0.569444444,0.5,0.472222222,0.472222222,0.486111111)


# Method 1.
# https://rpsychologist.com/d3/equivalence/
# http://daniellakens.blogspot.com/2016/12/tost-equivalence-testing-r-package.html
# install.packages("TOSTER")
library("TOSTER")
TOSTone(mean(x), .5, sd(x), length(x), -.5, .5, .05 , plot = TRUE, verbose = TRUE)
# +/- .5, that's a medium Cohen's d effect size.

# Equivalence Test Result:
# The equivalence test was significant, t(15) = 1.848, p = 0.0422, given equivalence bounds of -0.0228 and 0.0228 (on a raw scale) and an alpha of 0.05.
# 
# Null Hypothesis Test Result:
# The null hypothesis test was non-significant, t(15) = -0.152, p = 0.881, given an alpha of 0.05.
# 
# Based on the equivalence test and the null-hypothesis test combined, we can conclude that the observed effect is statistically not different from zero and statistically equivalent to zero.


# Method 2.
# https://rdrr.io/cran/equivUMP/#vignettes
# install.packages("equivUMP")
library("equivUMP")
equiv.test(x, mu=.5)

# One Sample equivalence test
# 
# data:  x
# t = -0.15226, df = 15, ncp = 4, p-value = 4.26e-05
# alternative hypothesis: equivalence
# null values:
#   lower upper
# [1,]  -Inf    -1
# [2,]     1   Inf
# sample estimates:
#   d 
# -0.03806567 