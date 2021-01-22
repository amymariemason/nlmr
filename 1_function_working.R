setwd("./Matt Arnold Code")
require(methods)
require(parallel)
require(metafor)
require(matrixStats)
require(survival)
source("nlmr_functions.r")
source("nlme_summ_aes MA.r")

#load test data & attach
dat = read.csv("NEUexp_SCZout_MRdata.csv",  header=TRUE, stringsAsFactors=FALSE)
attach(dat)

# add random mean data

dat$mean<-rexp(n=nrow(dat),rate=2)

# use function

frac_poly_summ_mr(bx=dat$BetaXG,by=dat$BetaYG, bxse=dat$seBetaXG, byse=dat$seBetaYG, xmean=dat$mean)

