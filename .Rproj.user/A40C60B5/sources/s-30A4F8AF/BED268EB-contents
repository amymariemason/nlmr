###
#Author: Amy Mason
# Purpose: generate simulatoion data from nlmr paper 
# Date: Jan 2019
##
1
# creating random underlying data
# We want: 
  # g - a binary gene
  # u - an unmeasure confounder
  # use these to create x (the exposure e.g. BMI) and y (the outcome e.g. blood pressure)
  # errorX is the error term on X
  # errorY is a random error term on Y
#In line with nlmr paper take:G~Bin(2,0.3), U~Uni(0,1), Ex~Exp(1), Ey~ N(0,1) 

data$u<-runif(1000,0,1)
data$errorX<-rexp(10000,1)
data$errorY<-rnorm(10000,0,1)
data$X<-2+0.25*data$g+data$u +data$errorX
beta1<-1.5
beta2<-0.5

# generate various Y with different exposure-outcome results
data$linear.Y<-beta1*data$X+0.8*data$u+data$errorY
data$quadratic.Y<-beta1*(data$X)^2+beta2*data$X+0.8*data$u+data$errorY
data$sqrt.Y<-beta1*sqrt(data$X)+0.8*data$u+data$errorY
data$log.Y<-beta1*log(data$X)+0.8*data$u+data$errorY
data$threshold.Y<-beta1*data$X+0.8*data$u+data$errorY