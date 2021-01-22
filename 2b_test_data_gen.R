#############################################################
# test data pre-summary


data<-create_data(10000)
# TEST 1: Linear data
test1<- fracpoly_mr(data$linear.Y, data$X, data$g, family="gaussian", q = 10,
                    d = 1, fig = T)
test2<- piecewise_mr(data$linear.Y, data$X, data$g, family="gaussian", fig = T)

f <- function(x) (3*(x-mean(data$X)))
test1$figure+stat_function(fun = f, colour = "green")+
  test2$figure+stat_function(fun = f, colour = "green")
# TEST 2: quadratic data

test1<-  fracpoly_mr(data$quadratic.Y, data$X, data$g, family="gaussian",
                     q = 10, d = "both", fig = T)

test2<-  piecewise_mr(data$quadratic.Y, data$X, data$g, family="gaussian",
                      fig = T)

xref=mean(data$X)
f <- function(x) (7*x^2 + 3*x - 7*xref^2-3*xref)
test1$figure+stat_function(fun = f, colour = "green")+
  test2$figure+stat_function(fun = f, colour = "green")


# 
fracpoly_mr(data$sqrt.Y, data$X, data$g, family="gaussian", q = 10, d = 1,
            fig=T)
fracpoly_mr(data$log.Y, data$X, data$g, family="gaussian", q = 10, d = 1, 
            fig = T)
fracpoly_mr(data$threshold.Y, data$X, data$g, family="gaussian", q = 10, d = 1,
            fig = T)

# Conclusion: all data is generating correctly and identified by the programs
# for linear and quadratic data


##########################################################################
# test data post summary

summ<-summary_function(data=data, gene = "g", exposure = "X",
                       outcome = "linear.Y", quant = 10)

head(test1$xcoef)


# these agree on the betaXG term; close on BetaXY term (why not identical?)

