###
#Author: Amy Mason
# Purpose: test quadratic simulated data
# Date: Jan 2019
##
setwd("C:/Users/Amy/Documents/Non-linear MR/Matt Arnold Code")
require(MASS)
require(methods)
require(parallel)
require(metafor)
require(ggplot2)
require(matrixStats)
require(survival)
require(ggplot2)
require(gridExtra)
source("nlmr_functions.r")
source("nlme_summ_aes.r")
# fixing random numbers for repetition of what is generated
set.seed(4743045)
# creating random underlying data


###########################################################
# test nlme_summ_aes on summerised data (1 set)

beta_set1<-1.5
beta_set2<-30

lotsofdata<- create_summary_data(Ytype = "quad", keep = TRUE, N = 10000, 
                                 beta1 = beta_set1, beta2 = beta_set2,
                                 quantiles = 100, confound = 0.5)
testdata<-lotsofdata$summary 
alldata<-lotsofdata$alldata

# plot the underlying data entire

# by genetic type
ggplot(data = alldata, aes(x = X, y = "quadratic.Y") )+
  geom_jitter(alpha = 0.3, aes(colour = as.factor(g)))

# underlying data check
alldata$xsq<-alldata$X^2
lm(data = alldata, quadratic.Y~X+xsq)


# non summerised fracpoly
keep1 <- fracpoly_mr(alldata$quadratic.Y, alldata$X, alldata$g, family = "gaussian",
            q = 10, d = 1, fig = T)


# summerised
keep2 <- frac_poly_summ_mr(bx = testdata$BetaXG, bxse = testdata$seBetaXG,
                        by = testdata$BetaYG, byse = testdata$seBetaYG,
                        xmean = testdata$meanX, family ="gaussian",
                        fig = TRUE, d = "both")
summary.frac_poly_mr(keep2)

f <- function(x) (beta_set1*x +
                    beta_set2*x^2 - beta_set1*3.65 -
                    beta_set2*3.65^2 )

plot1 <- keep1$figure+ stat_function(fun = f, colour = "green")+
  labs(title="fracpoly_mr")
plot2 <- keep2$figure+ stat_function(fun = f, colour = "green") +
  labs(title="fracpoly_summ_mr")
grid.arrange(plot1, plot2, ncol=2)

#######

# linear.Y = beta1*X
# quad Y = <-beta2*X^2+beta1*X
# sqrt.Y = beta1*sqrt(X)
# log.Y = beta1*log(X)
# threshold.Y<-ifelse(X>beta2,beta1*X,0)
