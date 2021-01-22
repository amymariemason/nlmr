###
#Author: Amy Mason
# Purpose: functions to generate simulation data from nlmr paper 
# Date: Jan 2019
##
if (!require("pacman")) install.packages("pacman")
pacman::p_load(MASS,methods,parallel, metafor, ggplot2, matrixStats, survival,
               patchwork)
if (!require("remotes")) install.packages("remotes")
if (!require("nlmr")) remotes::install_github("jrs95/nlmr")
library(nlmr)
source("../Matt Arnold Code/nlmr_functions.r")
source("../Matt Arnold Code/nlme_summ_aes.R")
# fixing random numbers for repetition of what is generated in test
#set.seed(4743045)
# creating random underlying data

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


create_data <- function(N, p1=1, p2=0, beta0=0, beta1 = 3, beta2 = 7, confound = 0.8){
  # generate G  
  data <- as.data.frame(rbinom(N,2,0.3))
  names(data) <- c("g")
  
  # generate Unknown confound
  data$u <- runif(N,0,1)
  
  # generate error terms
  data$errorX <- rexp(N,1)
  data$errorY <- rnorm(N,0,1)
  
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
  if(p1==p2){
    if(p1==0){
    data$fracpoly.Y <- beta0+beta1*log(data$X)+beta2*log(data$X)*log(data$X) +
                        confound*data$u + data$errorY
    }else{
      data$fracpoly.Y <- beta0+beta1*data$X^p1+beta2*log(data$X)*data$X^p1 +
        confound*data$u + data$errorY  
      }
    }else{
      if(p1==0){ 
        data$fracpoly.Y <- beta0+beta1*log(data$X) + beta2*data$X^p2 +
        confound*data$u + data$errorY
        
      }else if(p2==0){
        data$fracpoly.Y <-beta0+beta1*data$X^p1 + beta2*log(data$X) +
          confound*data$u + data$errorY
      }else{
    data$fracpoly.Y <- beta0+beta1*data$X^p1 + beta2*data$X^p2 +
                        confound*data$u + data$errorY
      }
    }
  return(data)
}

#############################################################
#  create quanta summary data function
#  some code recycled from Stephen Burgess
 
 
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
 
 
 
summary_function <- function(data, gene, exposure, outcome, quant = 10) { 

  ivf<- iv_free(y = data[,outcome], x = data[,exposure], g = data[,gene],
                covar = NULL, q = quant, family = "gaussian") 
  
  quantx0 <- ivf$x0q
  
  # this calculates the association for each quanta
  BetaYG   = NULL
  seBetaYG = NULL
  BetaXG  = NULL
  seBetaXG = NULL
  meanX = NULL

  for (j in 1:quant) {
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
# print(list(summary = head(output)))
 invisible(list(summary = output))
}
 


#############################################################
# generate summary data

create_summary_data<-function(Ytype = "linear", quantiles = 10, keep = FALSE,
                             ...) {
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
  } else if(Ytype == "fracpoly"){ 
    Ytype_name <- "fracpoly.Y"  
  } else {
    stop("model type not supported")

      }

  # create the data
  
  data<-create_data(...)
  summ<-summary_function(data, gene = "g", exposure = "X",
                         outcome = Ytype_name, quant = quantiles)
  summ_data<-summ$summary
    
  # keep entire set if keep variable set to TRUE
   if(keep == TRUE){
      data$quantiles <- summ$quantilesort
 #     print(paste(N, "individuals generated and summerized into ", quantiles,
 #                 " quantiles"))
 #     print("individual data kept")
 #     print(list(summary = head(summ_data), alldata = head(data)))
      invisible(list(summary = summ_data, alldata = data))
   }else{
#      print(paste(N, "individuals generated and summerized into ", quantiles,
#                  " quantiles"))
 #     print("individual data not kept")
#      print(list(summary = head(summ_data)))
      invisible(list(summary = summ_data))
  }
  }

