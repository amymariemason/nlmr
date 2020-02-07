setwd("U:/My Documents/Non-linear MR/Matt Arnold Code")
require(methods)
require(parallel)
require(metafor)
require(ggplot2)
require(matrixStats)
require(survival)
source("nlmr_functions.r")
source("nlme_summ_aes MA.r")
setwd("U:/My Documents/Code_review/Foley_2019/CodeReview1/")
#load test data & attach
dat = read.table("NEUexp_SCZout_MRdata.txt", sep=" ", header=TRUE, stringsAsFactors=FALSE)
attach(dat)

# add random mean data

dat$mean<-rexp(n=nrow(dat),rate=2)

# use function

frac_poly_summ_mr(bx=BetaXG,by=BetaYG, bxse=seBetaXG, byse=seBetaYG, xmean=dat$mean)

