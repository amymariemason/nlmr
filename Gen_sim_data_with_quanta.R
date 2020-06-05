###
#Author: Amy Mason
# Purpose: generate simulatoion data from nlmr paper 
# Date: Jan 2019
##
library(MASS)
# fixing random numbers for repetition of what is generated
set.seed(4743045)
# creating random underlying data

create_data <- function(N, beta1=1.5, beta2=0.5) {
# generate G  
data<-as.data.frame(rbinom(N,2,0.3))
names(data)<-c("g")

# generate U,
data$u<-runif(N,0,1)

# generate error terms
data$errorX<-rexp(N,1)
data$errorY<-rnorm(N,0,1)

# build X
data$X<-2+0.25*data$g+data$u +data$errorX

# generate various Y with different exposure-outcome results
data$linear.Y<-beta1*data$X+0.8*data$u+data$errorY
data$quadratic.Y<-beta1*(data$X)^2+beta2*data$X+0.8*data$u+data$errorY
data$sqrt.Y<-beta1*sqrt(data$X)+0.8*data$u+data$errorY
data$log.Y<-beta1*log(data$X)+0.8*data$u+data$errorY
data$threshold.Y<-ifelse(data$X>beta2,beta1*data$X,0)+0.8*data$u+data$errorY
return(data)
}

# function: create_data
# creates a dataset of individuals with a genetically influences outcome and exposure
# inputs
  # N number of individuals to generate
  # beta1, beta2: parameters for the 5 functions
# outputs 
  # g - a binary gene
  # u - an unmeasure confounder
  # use these to create x (the exposure e.g. BMI) and y (the outcome e.g. blood pressure)
  # errorX is the error term on X
  # errorY is a random error term on Y
#In line with nlmr paper take:G~Bin(2,0.3), U~Uni(0,1), Ex~Exp(1), Ey~ N(0,1) 

#############################################################
# test data pre-summary

# library(devtools)
# install_github("jrs95/nlmr")


library(nlmr)
data<-create_data(10000)
 fracpoly_mr(data$linear.Y, data$X, data$g, family="gaussian", q=10, d=1, fig=T)
 fracpoly_mr(data$quadratic.Y, data$X, data$g, family="gaussian", q=10, d=1, fig=T)
 fracpoly_mr(data$sqrt.Y, data$X, data$g, family="gaussian", q=10, d=1, fig=T)
 fracpoly_mr(data$log.Y, data$X, data$g, family="gaussian", q=10, d=1, fig=T)
 fracpoly_mr(data$threshold.Y, data$X, data$g, family="gaussian", q=10, d=1, fig=T)

# 
 
#############################################################
#  create quanta summary data function
summary_function <- function(data, gene, exposure, outcome) { # linear model of G->Y
  YG<-lm(data[,outcome]~data[,gene])
  BetaYG<-summary(YG)$coefficients[2,1]
  seBetaYG<-summary(YG)$coefficients[2,2]    
# linear model of G->X
XG<-lm(data[,exposure]~data[,gene])
BetaXG<-summary(XG)$coefficients[2,1]
seBetaXG<-summary(XG)$coefficients[2,2]
# mean of X
meanX<-mean(data[,exposure])
 output<-c(BetaXG, BetaYG, seBetaXG, seBetaYG, meanX)
	return(output)
}
# function: summary_function
# creates the needed quantile summary values from a generated data set
# inputs
  # gene: vector of presence/absense of gene
  # exposure: vector of exposure values
  # outcome: vector of outcome values
  # number of quanta
# outputs
  # BetaXG[j] :  association between X and G in quantile j
  # BetaYG[j] : association between Y and G in quantile j
  # seBetaXG[j] : s.e. for G-> X in quantile j
  # seBetaYG[j] :  standard error for G-> Y in quantile j
  # meanX[j]: average value of X in quantile j
 
 
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

#############################################################
# generate summary data

create_summary_data<-function(N1, Ytype="linear", keep=FALSE, N2=10000, beta1=1.5, beta2=0.5) {
  # create empty summary set
  df <- data.frame(matrix(ncol = 6, nrow = 0))
  names(df) <- c("Set", "BetaXG", "BetaYG", "seBetaXG", "seBetaYG", "meanX")
  
  # create Ytype_name to generate the appropiete function type
    if(Ytype=="linear"){
      Ytype_name<-"linear.Y"
    }
    if(Ytype=="quad"){
      Ytype_name<-"quadratic.Y"
    }
    if(Ytype=="sqrt"){
      Ytype_name<-"sqrt.Y"
    }
    if(Ytype=="log"){ 
      Ytype_name<-"log.Y"
    }
    if(Ytype=="threshold"){ 
      Ytype_name<-"threshold.Y"
    }
   if(!Ytype %in% c("linear", "quad", "sqrt", "log", "threshold")){
      stop("model type not supported")
   }
  
  # loop over N1: number of sets to make
  
  for(n in 1:N1) {
    #print(n)
    data<-create_data(N2, beta1, beta2)
    
     # sometimes this fails if the dataset generated is small and the linear models can't be fitted
    # error catch to simply regenerate data  
      newline <- try({
        cbind(n,t(summary_function(data, gene="g", exposure="X", outcome=Ytype_name)))
        }, silent=TRUE)
    # add error counter to force end of function if this is happening more than 10% of the time
      error_count=0
      give_up= FALSE
    #  
      while (class(newline)=="try-error" & give_up==FALSE){
        error_count = error_count+1
        if ((n>100 & error_count/n<0.5)|(n>100 & error_count/n<0.1)){give_up=TRUE}
        data<-create_data(N2, beta1, beta2)
        newline <- try({
          cbind(n,t(summary_function(data, gene="g", exposure="X", outcome="linear.Y")))
          }, silent=TRUE)
      }
      
      # if fails more than is reasonable
       if (give_up==TRUE){ 
        stop("Summarizing data in model failed; recommend changing data generation settings.")
        }
      
   if(class(newline)!="try-error"){
    names(newline) <- c("Set", "BetaXG", "BetaYG", "seBetaXG", "seBetaYG", "meanX")
    df<-rbind(df, newline )
  }
    
     # keep entire set if keep variable set to TRUE
  if(keep==TRUE){
    # append old set to end of this on subsequent rounds
    if(n==1){
      savedata<-data
    }
    else{
      savedata<-rbind(savedata,data)
    }
  } 
  # end of N1 loop      
  }
  #check sufficient summary sets are made
  if(nrow(df)<N1){stop(paste0("error: program terminated before ", N1, "datasets generated"))}
  
   names(df) <- c("Set", "BetaXG", "BetaYG", "seBetaXG", "seBetaYG", "meanX")
   # output data
   if(keep==TRUE){
      print(paste(N1, "datasets generated and summerized"))
      print(list(summary=head(df), alldata=head(savedata)))
     invisible(list(summary=df, alldata=savedata))
   }
    else{
      print(paste(N1, "datasets generated and summerized"))
      print(head(df))
      invisible(list(summary=df))
  }
  }



# function: create_summary_data
  
#output
    # data data.frame of summary statistics from generated sets of individual data
    # all data: copy of the individual level data, indexed by generation
  
break
###########################################################
# test nlme_summ_aes on summerised data

lotsofdata<- create_summary_data(10000, Ytype="quad", keep=TRUE, N2=50, beta1=1.5, beta2=0.5)
testdata<-lotsofdata$summary 

# plot the underlying data entire

ggplot(data=data$alldata, aes(x=X, y=linear.Y) )+ geom_jitter(alpha=0.3, aes(colour=as.factor(g)))+facet_wrap(facets=as.factor(data$alldata$g))

# plot the estimated effect of x on y, along the underlying model, via the frac_poly_summ_mr( code

keep<-frac_poly_summ_mr(bx=testdata$BetaXG,bxse=testdata$seBetaXG, by=testdata$BetaYG, byse=testdata$seBetaYG, xmean=testdata$meanX, family="gaussian",fig = TRUE, d="both")
summary.frac_poly_mr(keep)

f <- function(x) (0.5*x^2+1.5*x)
keep$figure+stat_function(fun=f, colour="green") 
 
# error tests
 
 by = testdata$BetaYG
 bx = testdata$BetaXG
 byse =testdata$seBetaYG
 bxse =bxse=testdata$seBetaXG
 xmean = testdata$meanX
 method="FE"
 d=1
 powers=c(0, -3, -2, -1.5, -1, -0.5, 1, 2)
 pd=0.05
 ci="model_se"
 nboot=100
 fig=FALSE
 family="binomial"
 offset=0
 pref_x="x"
 pref_y="y"
 ref=NA
 ci_type="overall"
 breaks=NULL
 ylim_lower = NA
 ylim_upper = NA
 xlim_lower = NA
 xlim_upper = NA