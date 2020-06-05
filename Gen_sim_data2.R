###
#Author: Amy Mason
# Purpose: generate simulation data from nlmr paper 
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

create_data <- function(N, beta1 = 1.5, beta2 = 0.5, confound = 0.8,
                        errorvar = 1) {
  # generate G  
  data <- as.data.frame(rbinom(N,2,0.3))
  names(data) <- c("g")
  
  # generate U,
  data$u <- runif(N,0,1)
  
  # generate error terms
  data$errorX <- rexp(N,errorvar)
  data$errorY <- rnorm(N,0,errorvar)
  
  # build X
  data$X<-2 + 0.25*data$g + data$u + data$errorX
  
  # generate various Y with different exposure-outcome results
  data$linear.Y <-      beta1*data$X + confound*data$u + data$errorY
  data$quadratic.Y <-   beta2*(data$X)^2 + beta1*data$X + confound*data$u +
                          data$errorY
  data$sqrt.Y <-        beta1*sqrt(data$X) + confound*data$u + data$errorY
  data$log.Y <-         beta1*log(data$X) + confound*data$u + data$errorY
  data$threshold.Y <-   ifelse(data$X > beta2, beta1*data$X, 0)+
                          confound*data$u + data$errorY
  return(data)
}

# function: create_data
# creates a dataset of individuals with a genetically influences outcome and exposure
# inputs
  # N number of individuals to generate
  # beta1, beta2: parameters for the 5 functions
# outputs 
  # g - a binary gene
  # u - an unmeasured confounder
  # use these to create x (the exposure e.g. BMI) and y (the outcome e.g. blood pressure)
  # errorX is the error term on X
  # errorY is a random error term on Y
#In line with nlmr paper take:G~Bin(2,0.3), U~Uni(0,1), Ex~Exp(1), Ey~ N(0,1) 

### POTENTIAL EXPANSION - add known covariate (linear?)

#############################################################
# test data pre-summary

library(nlmr)
data<-create_data(10000)
 fracpoly_mr(data$linear.Y, data$X, data$g, family="gaussian", q = 10, d = 1,
             fig = T)
 fracpoly_mr(data$quadratic.Y, data$X, data$g, family="gaussian", q = 10, d = 1,
             fig = T)
 fracpoly_mr(data$sqrt.Y, data$X, data$g, family="gaussian", q = 10, d = 1,
             fig=T)
 fracpoly_mr(data$log.Y, data$X, data$g, family="gaussian", q = 10, d = 1, 
             fig = T)
 fracpoly_mr(data$threshold.Y, data$X, data$g, family="gaussian", q = 10, d = 1,
             fig = T)

# 
 
#############################################################
#  create quanta summary data function
#  some code recycled from Stephen Burgess
summary_function <- function(data, gene, exposure, outcome, quant = 100) { 
  # linear model of G->Y
  YG<-lm(data[,outcome]~data[,gene])
  Y0<- data[,outcome] - YG$fit
  
  # this calculates the "IV-free exposure"
  qs = quantile(Y0, prob=seq(0, 1-1/quant, by=1/quant))
  # this divides into strata based on IV-free exposure
  quantx0 = as.numeric(lapply(Y0, function(x) { return(sum(x>qs)) }))

  # this calculates the association for each quanta
  BetaYG   = NULL
  seBetaYG = NULL
  BetaXG  = NULL
  seBetaXG = NULL
  meanX = NULL
  
  for (j in 1:length(qs)) {
    BetaYG[j]   = summary(lm(data[quantx0 == j, outcome] ~ 
                               data[quantx0 == j, gene]))$coef[2]
    seBetaYG[j] = summary(lm(data[quantx0 == j,outcome] ~
                               data[quantx0 == j, gene]))$coef[2,2]
    BetaXG[j]   = summary(lm(data[quantx0 == j, exposure] ~ 
                               data[quantx0 == j, gene]))$coef[2]
    seBetaXG[j] = summary(lm(data[quantx0 == j, exposure] ~
                               data[quantx0 == j, gene]))$coef[2,2]
    meanX[j] = mean(data[quantx0 == j, exposure])
  }
  

 output <- data.frame(BetaXG, BetaYG, seBetaXG, seBetaYG, meanX)
 print(list(summary = head(output)))
 invisible(list(summary = output))
}
 
# function: summary_function
# creates the needed quantile summary values for an nlmr from a individual level data set
# inputs
  # gene: vector of presence/absense of gene
  # exposure: vector of exposure values
  # outcome: vector of outcome values
  # quant: number of quantiles
# outputs: data.frame containing
  # BetaXG[j] :  association between X and G in quantile j
  # BetaYG[j] : association between Y and G in quantile j
  # seBetaXG[j] : s.e. for G-> X in quantile j
  # seBetaYG[j] :  standard error for G-> Y in quantile j
  # meanX[j]: average value of X in quantile j
 
 
# NOTE: the stratification is done on residual phenotype 
# (residual phenotype  = individual's phenotype - centred genetic contribution 
# to phenotype from included genetic variants)
# this means we compare individuals in the population who would have similar 
# phenotype value IF they had the same genetic code
# see https://www.bmj.com/content/364/bmj.l1042 for details of how this was 
# done with BMI

 


#############################################################
# generate summary data

create_summary_data<-function(Ytype = "linear", quantiles = 100, keep = FALSE,
                              N = 10000, beta1 = 1.5, beta2 = 0.5, 
                              confound = 0.8, errorvar = 1) {
  # create empty summary set
  df <- data.frame(matrix(ncol = 6, nrow = 0))
  names(df) <- c("Set", "BetaXG", "BetaYG", "seBetaXG", "seBetaYG", "meanX")
  
  # create Ytype_name to generate the appropiete function type
  if(Ytype == "linear"){
    Ytype_name <- "linear.Y"
  } else if (Ytype == "quad"){
    Ytype_name <- "quadratic.Y"
  } else if(Ytype == "sqrt"){
    Ytype_name <- "sqrt.Y"
  } else if(Ytype == "log"){ 
    Ytype_name <- "log.Y"
  } else if(Ytype == "threshold"){ 
    Ytype_name <- "threshold.Y"
  } else {
    stop("model type not supported")
  }

  # create the data
  
  data<-create_data(N, beta1, beta2, confound)
  summ<-summary_function(data, gene = "g", exposure = "X",
                         outcome = Ytype_name, quant = quantiles)
  summ_data<-summ$summary
    
  # keep entire set if keep variable set to TRUE
   if(keep == TRUE){
      data$quantiles <- summ$quantilesort
      print(paste(N, "individuals generated and summerized into ", quantiles,
                  " quantiles"))
      print("individual data kept")
      print(list(summary = head(summ_data), alldata = head(data)))
      invisible(list(summary = summ_data, alldata = data))
   }else{
      print(paste(N, "individuals generated and summerized into ", quantiles,
                  " quantiles"))
      print("individual data not kept")
      print(list(summary = head(summ_data)))
      invisible(list(summary = summ_data))
  }
  }



# function: create_summary_data
  
#output
    # data data.frame of summary statistics from generated sets of individual data
    # all data: copy of the individual level data
  

###########################################################
# test nlme_summ_aes on summerised data

lotsofdata<- create_summary_data(Ytype = "quad", keep = TRUE, N = 100000, 
                                 beta1 = 1.5, beta2 = 0.5, quantiles = 100,
                                 confound = 0)
testdata<-lotsofdata$summary 
alldata<-lotsofdata$alldata
# recentre meanX
testdata$meanX<-testdata$meanX -mean(alldata$X)
# plot the underlying data entire

# by genetic type
ggplot(data = alldata, aes(x = X, y = "linear.Y") )+
  geom_jitter(alpha = 0.3, aes(colour = as.factor("g")))+
  facet_wrap(facets = as.factor(alldata$g))


# by residual strata
precision<-15
roundUp <- function(x) ceiling(x / precision) * precision
alldata$decile<-ifelse(alldata$quantile>0,roundUp(alldata$quantiles), precision)
ggplot(data=alldata, aes(x=X, y=quadratic.Y) )+ geom_jitter(alpha=0.3, aes(colour=as.factor(decile)))+facet_wrap(facets=as.factor(alldata$decile))


ggplot(data=alldata, aes(x=X, y=quadratic.Y) )+ geom_jitter(alpha=0.3, aes(colour=as.factor(decile)))+facet_wrap(facets=as.factor(alldata$g))

# plot the estimated effect of x on y, along the underlying model, via the frac_poly_summ_mr( code

keep<-frac_poly_summ_mr(bx=testdata$BetaXG,bxse=testdata$seBetaXG, by=testdata$BetaYG, byse=testdata$seBetaYG, xmean=testdata$meanX, family="gaussian",fig = TRUE, d="both")
summary.frac_poly_mr(keep)

f <- function(x) (0.5*x^2+1.5*x)
keep$figure+stat_function(fun=f, colour="green") 
 
#########################################################
# STEVE HERE: just test
lotsofdata<- create_summary_data(Ytype = "quad", keep = TRUE, N = 10000,
                                 beta1 = 1.5, beta2 = 7, quantiles = 20,
                                 confound = 0, errorvar = 0)
testdata <- lotsofdata$summary 
alldata <- lotsofdata$alldata
testdata$meanX <- testdata$meanX

# underlying data check
alldata$sqX <- alldata$X^2
lm(data = alldata, quadratic.Y~X + sqX)

# non summerised fracpoly
fracpoly_mr(alldata$quadratic.Y, alldata$X, alldata$g, family = "gaussian",
            q = 10, d = 1, fig = T)


# summerised
keep<-frac_poly_summ_mr(bx = testdata$BetaXG, bxse = testdata$seBetaXG,
                        by = testdata$BetaYG, byse = testdata$seBetaYG,
                        xmean = testdata$meanX, family ="gaussian",
                        fig = TRUE, d ="both")
summary.frac_poly_mr(keep)

f <- function(x) (7*x^2 + 1.5*x)
keep$figure + stat_function(fun = f, colour = "green")

# linear.Y = beta1*X
# quad Y = <-beta2*X^2+beta1*X
# sqrt.Y = beta1*sqrt(X)
# log.Y = beta1*log(X)
# threshold.Y<-ifelse(X>beta2,beta1*X,0)
