###
#Author: Amy Mason
# Purpose: test linear simulated data
# Date: Jan 2019
##
setwd("../Matt Arnold Code")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(MASS,methods,parallel, metafor, ggplot2, matrixStats, survival,
               patchwork)
if (!require("remotes")) install.packages("remotes")
if (!require("nlmr")) remotes::install_github("jrs95/nlmr")
library(nlmr)
source("nlmr_functions.r")
source("nlme_summ_aes.R")
source("../2_Gen_sim_data.R")
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
            q = 10, d = "both", fig = F)


# summerised
keep2 <- frac_poly_summ_mr(bx = testdata$BetaXG, bxse = testdata$seBetaXG,
                        by = testdata$BetaYG, byse = testdata$seBetaYG,
                        xmean = testdata$meanX, family ="gaussian",
                        fig = F, d ="both")
summary.frac_poly_mr(keep2)

f <- function(x) (beta_set*(x)- beta_set*mean(alldata$X))


plot1 <- keep1$figure+ stat_function(fun = f, colour = "green") +
  ggtitle("fracpoly_mr")
plot2 <- keep2$figure+ stat_function(fun = f, colour = "green")+
  ggtitle("fracpoly_summ_mr")

plot1+plot2

