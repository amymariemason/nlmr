# will generate the degree 1 outcome of nlmr frac poly testing for specific setting
# author: Amy Mason
# input: number of repeats, par1, par2, beta1, beta2
# input: 




setwd("/rds/project/jmmh2/rds-jmmh2-projects/zz_mr/nlmr/Non-linear MR/Matt Arnold Code")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(MASS,methods,parallel, metafor, ggplot2, matrixStats, survival,
               patchwork)
if (!require("remotes")) install.packages("remotes")
if (!require("nlmr")) remotes::install_github("jrs95/nlmr")
library(nlmr)
source("../nlme_summ_aes_AM.R")
# source("nlme_summ_aes.r")
source("../2_Gen_sim_data.R")
source("../fit_power_AM.R")
source("../4_assess_output.R")

# slurm code to run 
# sbatch --export=basePath="a" 5_gen_save_output.R

par1<- as.numeric(Sys.getenv('par1'))
beta1<-as.numeric(Sys.getenv('beta1'))


# add parellel calculations

library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
set.seed(4743045)

# create the data
all<-  foreach(i=1:500,
          .inorder=FALSE,
          .packages = c("nlmr","metafor")) %dopar% {
            source("../4_assess_output.R")  
    checkEst(x=i,deg=1, par1=par1,par2=0,beta0=0,beta1=beta1,beta2=0)
    
          }
stopCluster(cl)
all2<-as.data.frame(do.call(rbind,all))
test1_all<-as.data.frame(do.call(rbind,all2$Test1_all))
test1_summ<-as.data.frame(do.call(rbind,all2$Test1_summ))

# summerise test result in single line  (d=1)

Answer_deg1<-c("alldata",nrow(test1_all),test1_all[1, c("p1",  "beta1")],
               
               
               mean(test1_all$Beta_1_est), sd(test1_all$Beta_1_est), mean(test1_all$Beta_1_se),
               mean(test1_all$Beta_1_cov), 
               mean(test1_all$powers_correct), mean(test1_all$powers_in_set)
)

names(Answer_deg1)<-c("type","simulations", "p1", "Beta1", 
                      "est_Beta1_mean", "est_Beta1_sd", "est_Beta1_mean_SE",
                      "cov_beta1", 
                      "fit_correct", "fit_set"
)

Answer_deg1_summ<-c("summ", nrow(test1_summ),test1_summ[1, c("p1",  "beta1")],
                    mean(test1_summ$Beta_1_est), sd(test1_summ$Beta_1_est), mean(test1_summ$Beta_1_se),
                    mean(test1_summ$Beta_1_cov), 
                    mean(test1_summ$powers_correct), mean(test1_summ$powers_in_set)
)

names(Answer_deg1_summ)<-c("type","simulations", "p1", "Beta1", 
                           "est_Beta1_mean", "est_Beta1_sd", "est_Beta1_mean_SE",
                           "cov_beta1", 
                           "fit_correct", "fit_set")


# look at the difference in p-test 

test1_diff<- test1_all[,c("fp_d1_d2", "fp", "quad","Q","cochQ","trend" )] - 
  test1_summ[,c("fp_d1_d2", "fp", "quad","Q","cochQ","trend" )]

max_diffs<-apply(as.matrix(abs(test1_diff)),2,max)
names(max_diffs)<- c("max_fp_d1_d2_diff", "max_fp_diff", "max_quad_diff",
                     "max_Q_diff","max_cochQ_diff","max_trend_diff")
mean_diffs<-apply(as.matrix(test1_diff),2,mean)
names(mean_diffs)<- c("mean_fp_d1_d2_diff", "mean_fp_diff", "mean_quad_diff",
                      "mean_Q_diff","mean_cochQ_diff","mean_trend_diff")
se_diffs<-apply(as.matrix(test1_diff),2,sd)
names(se_diffs)<- c("se_fp_d1_d2_diff", "se_fp_diff", "se_quad_diff",
                    "se_Q_diff","se_cochQ_diff","se_trend_diff")

# how often do these differences create a change in significance at the 5% sig level
test1_sig<-mapply(test1_all[,c("fp_d1_d2", "fp", "quad","Q","cochQ","trend" )], FUN=function(x)ifelse(x<0.05,1,0))
test1_summ_sig <-mapply(test1_summ[,c("fp_d1_d2", "fp", "quad","Q","cochQ","trend" )], FUN=function(x)ifelse(x<0.05,1,0))
sig_diff<-abs(test1_sig-test1_summ_sig)

sig_changes<- apply(sig_diff,2, mean)
names(sig_changes)<- c("change_fp_d1_d2_5%", "change_fp_5%", "change_quad_5%",
                       "change_Q_5%","change_cochQ_5%","change_trend_5%")

Ans_deg1<-rbind(Answer_deg1, Answer_deg1_summ)
Ans_deg1_b<-rbind(c(max_diffs, mean_diffs, se_diffs, sig_changes), 
                  c(max_diffs, mean_diffs, se_diffs, sig_changes))
Ans1<-cbind(Ans_deg1,Ans_deg1_b)
Ans1<-as.data.frame(Ans1)
Ans1$date<-as.character(Sys.time())
Ans1$date2<-Sys.time()
Ans1<-as.matrix(Ans1)




# write data out to file

#check if file exists - if not makes headers 


  names_out<- c("type","simulations", "p1", "Beta1", 
                "est_Beta1_mean", "est_Beta1_sd", "est_Beta1_mean_SE",
                "cov_beta1", 
                "fit_correct", "fit_set",
                "max_fp_d1_d2_diff", "max_fp_diff", "max_quad_diff",
                "max_Q_diff","max_cochQ_diff","max_trend_diff", 
                "mean_fp_d1_d2_diff", "mean_fp_diff", "mean_quad_diff",
                "mean_Q_diff","mean_cochQ_diff","mean_trend_diff",
                "se_fp_d1_d2_diff", "se_fp_diff", "se_quad_diff",
                "se_Q_diff","se_cochQ_diff","se_trend_diff",
                "change_fp_d1_d2_5%", "change_fp_5%", "change_quad_5%",
                "change_Q_5%","change_cochQ_5%","change_trend_5%", "date", "date2")
  
if (!file.exists("../Table_output/table1.csv")) {  
write.table(x=t(names_out), file="../Table_output/table1.csv", col.names=FALSE, row.names=FALSE, sep=",")
  
}

write.table(x=Ans1, file="../Table_output/table1.csv", append=TRUE, col.names = FALSE, row.names=FALSE, sep=",")

