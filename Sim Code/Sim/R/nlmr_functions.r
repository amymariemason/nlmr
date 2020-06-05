egfr <- function(creat, age, sex) {

	creat_over_k <- creat / ifelse(toupper(sex)=="MALE",0.9,0.7)
	
	min_part <- ifelse(creat_over_k < 1, creat_over_k, 1)
	max_part <- ifelse(creat_over_k > 1, creat_over_k, 1)
	
	res <- 141 * (min_part^ifelse(toupper(sex)=="MALE",-0.411,-0.329)) * (max_part^(-1.209)) * (0.993^age) * ifelse(toupper(sex)=="MALE",1,1.018)
	
	return(res)

}



# Liam: need to edit this function so that the coxph function runs the correct 
# Cox regression for the MVP scenario
grs_outcome <- function(q_i, q_var, indata, outcome, timevar) {

	data_q <- indata[indata[[q_var]] == q_i,]
	# standardise names
	data_q$timevar <- data_q[[timevar]]
	data_q$outcome <- data_q[[outcome]]
	
	cox_out <- coxph(Surv(time=timevar, event=outcome) ~ egfr_grs_centre + ages + as.factor(sex) + offset(wt_pre) + cluster(id) + strata(centreEH), data=data_q)
	
	grs_row <- grep("egfr_grs_centre", row.names(summary(cox_out)$coef))
	
	return(summary(cox_out)$coef[grs_row,c(1,4)])
	

}
# function to generate GRS - exposure betas and SEs 
# SK suggested fixed centreEH effect
grs_exposure <- function(q_i, q_var, indata) {
	data_q <- indata[indata[[q_var]] == q_i,]
	tmp <- lm(egfr ~ egfr_grs_centre + ages + ages^2 + as.factor(sex) + as.factor(centreEH), data=data_q)
	grs_row <- grep("egfr_grs_centre", row.names(summary(tmp)$coef))
	return(summary(tmp)$coef[grs_row,1:2])
}
# function to combine the GRS - exposure betas across strata to 
# get at a single "GRS -> eGFR" effect across the range of eGFR
est_overall_grs_exposure <- function(grs_exposure_mean_se) {
	return(sum(grs_exposure_mean_se[1,] * grs_exposure_mean_se[2,]^(-2)) / sum(grs_exposure_mean_se[2,]^(-2)))
}
# function to compute LACE estimates from the GRS - outcome associatons
# and the overall GRS - exposure association
lace_estimates <- function(grs_outcomes, grs_xcoef) {
	return(grs_outcomes / grs_xcoef)
}



## plot output 
## the complexity seems to be in making the piecewise linear sections join up correctly
# there's an arbitrariness in where the intercept of the overall curve is
# this means we can set the intercept in the first subset to be 0 
# and then shift everything later
redraw <- function(mean_se_mat) {

	out <- mean_se_mat
	out[1,] <- out[1,] + rnorm(dim(mean_se_mat)[2], 0, out[2,])
	return(out)

}	
calculate_piecewise_shape <- function(indata, q_var, xref, grs_outcome, grs_exposure) {

	q_index <- sort(unique(indata[[q_var]]))
	
	q_min <- sapply(q_index, function(z) min(indata$egfr[indata[[q_var]] == z]))
	q_max <- sapply(q_index, function(z) max(indata$egfr[indata[[q_var]] == z]))
	xmeans <- sapply(q_index, function(z) mean(indata$egfr[indata[[q_var]] == z]))
	
	grs_xcoef <- est_overall_grs_exposure(grs_exposure) 
	
	lace_e <- lace_estimates(grs_outcome, grs_xcoef)[1,]
	lace_se <- lace_estimates(grs_outcome, grs_xcoef)[2,]
	
	intercepts <- rep(0,length(q_index))
	for (q_i in q_index[-1]) {
		intercepts[q_i] <- intercepts[q_i - 1] + lace_e[q_i-1]*q_max[q_i-1] - lace_e[q_i]*q_max[q_i-1]
	}
	
	# calculate which group xref is in
	xref_q <- max(which(q_min < xref))

	# evaluate the value at the reference point
	xref_y <- intercepts[xref_q] + xref*lace_e[xref_q]

	# subtract this factor from all the intercepts
	intercepts <- intercepts - xref_y
	
	# calculate Y positions of the means
	means_y <- intercepts + lace_e*xmeans
	
	return(list(xmeans = xmeans, 
				grs_outcome = grs_outcome,
				grs_exposure = grs_exposure, 
				means_y = means_y, 
				intercepts = intercepts, 
				lace_e = lace_e, 
				lace_se = lace_se, 
				q_var=q_var,
				xref=xref,
				q_index = q_index, 
				q_min = q_min, 
				q_max = q_max))
	
}
plot_piecewise_shape <- function(shape_obj, nboot=5000, ncores=4) {

	
	# add a reference line
	abline(h=0, lty="dashed", col="grey80")
	
	# make all of the shape_obj values directly available in this function
	for (nam in names(shape_obj)) {
		assign(x=nam, value=shape_obj[[nam]])
	}
	
	# now plot the lines (do the first one specially)
	lines(x=c(q_min[1], q_max[1]), y=intercepts[1] + lace_e[1]*c(q_min[1], q_max[1]), lend=2)
	for (q_i in q_index[-1]) {
		lines(x=c(q_max[max(1,q_i-1)], q_max[q_i]), y=intercepts[q_i] + lace_e[q_i]*c(q_max[max(1,q_i-1)], q_max[q_i]), lend=2)
	}
	
	# bootstrap the whole piecewise shape calculation
	# technical steps required to make the data and functions 
	# available to worker nodes
	require(parallel)
	if (!exists("cl")) {
		cl <- makeCluster(ncores)
		clusterExport(cl, "calculate_piecewise_shape")
		clusterExport(cl, "redraw")
		clusterExport(cl, "est_overall_grs_exposure")
		clusterExport(cl, "lace_estimates")
		clusterExport(cl, "epic")
		clusterExport(cl, "shape_obj", envir=environment())
		
	}
	# run the bootstrap
	boot_reps <- parLapply(cl,1:nboot, function(i) { calculate_piecewise_shape(indata=epic, q_var=shape_obj$q_var, xref=shape_obj$xref, grs_outcome=redraw(shape_obj$grs_outcome), grs_exposure=shape_obj$grs_exposure)})
	# shutdown
	stopCluster(cl)
	rm(cl)
			
	y_means_boot <- (sapply(1:length(boot_reps), function(z) boot_reps[[z]]$means_y))
	cis <- apply(y_means_boot,1,quantile,probs=c(0.025,0.975))						

	for (q_i in q_index) {
		lines(x=rep(xmeans[q_i],2), y=cis[,q_i], lend=2, col="grey70")
	}

	# plot each of the means 
	points(xmeans, means_y, pch=16)
	

}		

## also use SB code to generate FP curves using the "summary" data
plot_fp_curve <- function(fit) {

	par(lend=2)
	lines(x=par("usr")[1:2], y=rep(0,2), col="grey80", lty="dashed")

	xpts_neworder <- order(fit$figure$data$x)
	xpts <- fit$figure$data$x[xpts_neworder]
	lci <- log(fit$figure$data$lci[xpts_neworder])
	uci <- log(fit$figure$data$uci[xpts_neworder])
	ypts <- log(fit$figure$data$yest[xpts_neworder])
	
	keep <- seq(from=1, to=length(xpts), by=2)

	points(xpts[keep], lci[keep], type="l", col="grey70")
	points(xpts[keep], uci[keep], type="l", col="grey70")
	points(xpts[keep], ypts[keep], type="l")


}					
