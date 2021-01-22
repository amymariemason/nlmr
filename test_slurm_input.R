base.path = Sys.getenv('basePath')
practise<-as.matrix( c( base.path, "test")) 
write.table(practise, "/rds/project/jmmh2/rds-jmmh2-projects/zz_mr/nlmr/Non-linear MR/Table_output/table_test.csv", append = TRUE)