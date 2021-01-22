install.packages(c("metafor", "ggplot2")) # only needs to be run once

library(metafor); library(ggplot2)
setwd("C:/Users/Amy/Documents/Non-linear MR/Matt Arnold Code") # replace with your project folder
source("nlme_summ_aes.R")     # this is the non-linear function

####
# generating the stratum specific summary estimates

pheno0  = pheno-lm(pheno~grs+sex+centre)$fit
          # add/remove covariates to taste
          # this calculates the "IV-free exposure"
quant = 100
          # we divided into 100 quantiles - can change
qs = quantile(pheno0, prob=seq(0, 1-1/quant, by=1/quant))

quantx0 = as.numeric(lapply(pheno0, function(x) { return(sum(x>qs)) }))
          # this divides into strata based on IV-free exposure

by   = NULL
byse = NULL
bx   = NULL
bxse = NULL
xmean = NULL

for (j in 1:length(qs)) {
 by[j]   = summary(coxph(outcome.surv[quantx0==j]~grs[quantx0==j]+sex[quantx0==j]))$coef[1,1]
 byse[j] = summary(coxph(outcome.surv[quantx0==j]~grs[quantx0==j]+sex[quantx0==j]))$coef[1,3]
 bx[j]   = summary(lm(pheno[quantx0==j]~grs[quantx0==j]+sex[quantx0==j]))$coef[2]
 bxse[j] = summary(lm(pheno[quantx0==j]~grs[quantx0==j]+sex[quantx0==j]))$coef[2,2]
 xmean[j] = mean(pheno[quantx0==j])
 }
 # this generates the stratum specific estimates
 # can replace coxph with glm or whatever
 # again, add covariates to taste

######################
# test dataset

data = read.csv("hunt_100_unweighted_acm.csv", header=TRUE)
 # this is a dataset we prepared previously with BMI as exposure

by   = data$beta_acm  # genetic associations with outcome per quantile
byse = data$se_acm    # standard errors of genetic associations with outcome
bx   = data$beta_bmi  # genetic associations with exposure per quantile
bxse = data$se_bmi    # standard errors of genetic associations with exposure
xmean = data$mean_bmi # mean exposure in each quantile

xmean1 = xmean-18  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
                   # if your exposure values start near zero, this isn't needed

# main figure

library(metafor); library(ggplot2)
frac_poly_summ_mr(by, bx, byse, bxse, xmean1, pd=0.5, ref=7, d=2, offset=18, xlim_upper=45, ylim_lower=NA, ylim_upper=10, fig=TRUE, pref_x="Body mass index", pref_y="Odds ratio for all-cause mortality", breaks=c(1,2.5,5,7.5,10))
 # lots of options here (see the help file ?frac_poly_summ_mr)
 # offset = 18 as we reduced BMI by 18
 # ref = 7 means that BMI = 25 is the reference category (25=7+18)
 # xlim and ylim control the axis limits
 # breaks is the labels on the y-axis

