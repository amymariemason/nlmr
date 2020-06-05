###
#Author: Amy Mason
# Purpose: test linear simulated data
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
source("nlmr_functions.r")
source("nlme_summ_aes MA.r")
# fixing random numbers for repetition of what is generated
set.seed(4743045)
# creating random underlying data


###########################################################
# test nlme_summ_aes on summerised data (1 set)

beta_set<-7

lotsofdata<- create_summary_data(Ytype = "linear", keep = TRUE, N = 1000, 
                                 beta1 = beta_set, beta2 = 0.5, quantiles = 30,
                                 confound = 0.5)
testdata<-lotsofdata$summary 
alldata<-lotsofdata$alldata

# plot the underlying data entire

# by genetic type
ggplot(data = alldata, aes(x = X, y = "linear.Y") )+
  geom_jitter(alpha = 0.3, aes(colour = as.factor(g)))

# underlying data check
lm(data = alldata, linear.Y~X)


# non summerised fracpoly
keep1 <- fracpoly_mr(alldata$linear.Y, alldata$X, alldata$g, family = "gaussian",
            q = 10, d = 1, fig = T)


# summerised
keep2 <- frac_poly_summ_mr(bx = testdata$BetaXG, bxse = testdata$seBetaXG,
                        by = testdata$BetaYG, byse = testdata$seBetaYG,
                        xmean = testdata$meanX, family ="gaussian",
                        fig = TRUE, d ="both")
summary.frac_poly_mr(keep2)

f <- function(x) (beta_set*(x- mean(alldata$X)))

require(gridExtra)
plot1 <- keep1$figure+ stat_function(fun = f, colour = "green") +
  title(main="fracpoly_mr")
plot2 <- keep2$figure+ stat_function(fun = f, colour = "green")+
  title(main="fracpoly_summ_mr")
grid.arrange(plot1, plot2, ncol=2)




#############################################

# linear.Y = beta1*X
# quad Y = <-beta2*X^2+beta1*X
# sqrt.Y = beta1*sqrt(X)
# log.Y = beta1*log(X)
# threshold.Y<-ifelse(X>beta2,beta1*X,0)
