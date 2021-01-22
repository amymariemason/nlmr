# fit poly

fracpoly_best_powers<-function (coef, coef_se, xmean, d = 1, pd = 0.05, method = "FE",
                                powers =  c(0, -3, -2, -1.5, -1, -0.5, 1, 2, 3)) 
{
  likelihood_d1 <- NULL
  for (p1 in powers) {
    if (p1 == -1) {
      x1 <- xmean^p1
    } else {
      x1 <- (p1 + 1) * xmean^p1
    }
    fp_mod <- try(rma(coef ~ -1 + x1, vi = (coef_se)^2, 
                      method = method), silent = TRUE)
    if (is(fp_mod, "try-error") == T) {
      likelihood_d1 <- c(likelihood_d1, NA)
    }  else {
      if (p1 == 0) {
        fp1 <- fp_mod
        p_ML <- p1
      }   else {
        if (fp_mod$fit.stats[1, 1] >= suppressWarnings(max(likelihood_d1, 
                                                           na.rm = T))) {
          fp1 <- fp_mod
          p_ML <- p1
        }
      }
      likelihood_d1 <- c(likelihood_d1, fp_mod$fit.stats[1, 
                                                         1])
    }
  }
  maxlik_d1 <- max(likelihood_d1, na.rm = T)
  fp_p <- 1 - pchisq(((-2 * likelihood_d1[1]) - (-2 * maxlik_d1)), 
                     df = 1)
  powers1 <- powers
  powers2 <- powers
  likelihood_d2 <- NULL
  for (p11 in powers1) {
    if (p11 == -1) {
      x1 <- xmean^p11
    }
    else {
      x1 <- (p11 + 1) * xmean^p11
    }
    for (p21 in powers2) {
      if (p11 == p21) {
        if (p21 == -1) {
          x2 <- 2 * (xmean^p21) * log(xmean)
        }
        else {
          x2 <- ((p21 + 1) * (xmean^p21) * log(xmean) + 
                   xmean^p21)
        }
      }
      else {
        if (p21 == -1) {
          x2 <- xmean^p21
        }
        else {
          x2 <- (p21 + 1) * xmean^p21
        }
      }
      fp_mod <- try(rma(coef ~ -1 + x1 + x2, vi = (coef_se)^2, 
                        method = method), silent = TRUE)
      if (is(fp_mod, "try-error") == T) {
        likelihood_d2 <- c(likelihood_d2, NA)
      }
      else {
        if (p11 == 0 & p21 == 0) {
          fp2 <- fp_mod
          p1_ML <- p11
          p2_ML <- p21
        }
        else {
          if (fp_mod$fit.stats[1, 1] >= suppressWarnings(max(likelihood_d2, 
                                                             na.rm = T))) {
            fp2 <- fp_mod
            p1_ML <- p11
            p2_ML <- p21
          }
        }
        likelihood_d2 <- c(likelihood_d2, fp_mod$fit.stats[1, 
                                                           1])
      }
    }
    powers2 <- powers2[-1]
  }
  maxlik_d2 <- max(likelihood_d2, na.rm = T)
  fp_d12_p <- 1 - pchisq(((-2 * maxlik_d1) - (-2 * maxlik_d2)), 
                         df = 2)
  if (d == "both") {
    if (fp_d12_p > pd) {
      d <- 1
    }
    else {
      d <- 2
    }
  }
  if (d == 1) {
    model <- fp1
    if (length(model$b) != 1) 
      stop("incorrect number of parameters for best fitting fractional polynomial of degree 1")
  }
  if (d == 2) {
    model <- fp2
    if (length(model$b) != 2) 
      stop("incorrect number of parameters for best fitting fractional polynomial of degree 2")
  }
  results <- list(model = model, p_ML = p_ML, p1_ML = p1_ML, 
                  p2_ML = p2_ML, fp_p = fp_p, fp_d12_p = fp_d12_p, d = d, loglik= model$fit.stats[1,1])
  return(results)
}

#######################################################################
# 

fit_poly<- function (y, x, g, covar = NULL, family = "gaussian", q = 10, 
          xpos = "mean", method = "FE", d = 1, pd = 0.05, ci = "model_se", 
          nboot = 100, fig = F, ref = mean(x), pref_x = "x", pref_x_ref = "x", 
          pref_y = "y", ci_type = "overall", ci_quantiles = 10, breaks = NULL,
          powers =  c(0, -3, -2, -1.5, -1, -0.5, 1, 2, 3) ) 
{
  if (!(is.vector(y) & is.vector(x) & is.vector(g))) 
    stop("either the outcome, exposure or instrument is not a vector")
  if (!is.null(covar)) {
    if (!is.data.frame(covar)) 
      stop("covar has to be a data.frame")
  }
  if (!((is.numeric(y) | is.integer(y)) & (is.numeric(x) | 
                                           is.integer(x)) & (is.numeric(g) | is.integer(g)))) 
    stop("either the outcome, exposure or instrument is not numeric")
  if (any(x <= 1)) 
    stop("fractional polynomial models require the exposure to be >>1")
  if (length(y) <= 1) 
    stop("the outcome is less than or equal to a single value")
  if (!(length(y) == length(x) & length(y) == length(g)) | 
      (if (!is.null(covar)) {
        (nrow(covar) != length(y))
      }
      else {
        FALSE
      })) 
    stop("the number of observations for the outcome, exposure, instrument and covariates are not all the same")
  if (any(is.na(y)) | any(is.na(x)) | any(is.na(g)) | (if (!is.null(covar)) {
    any(is.na(covar))
  }
  else {
    FALSE
  })) 
    stop("there are missing values in either the outcome, exposure, instrument or covariates")
  if (!(family == "gaussian" | family == "binomial")) 
    stop("family has to be equal to either \"gaussian\" or \"binomial\"")
  if (family == "binomial") {
    if (any(!(y == 1 | y == 0))) 
      stop("y has to be 0 or 1 if family is equal to \"binomial\"")
  }
  if ((length(y)/10) < q) 
    stop("the quantiles should contain at least 10 observations")
  if (!(xpos == "mean" | (xpos > 0 & xpos < 1))) 
    stop("the position used to relate x to the localised average causal effect")
  if (!(d == 1 | d == 2 | d == "both")) 
    stop("the degree has to be equal to 1, 2 or \"both\"")
  if (!(ci == "model_se" | ci == "bootstrap_se" | ci == "bootstrap_per")) 
    stop("the confidence intervals must be one of \"model_se\", \"bootstrap_se\" and \"bootstrap_per\"")
  if (!is.null(covar)) {
    covar <- model.matrix(as.formula(~.), data = covar)[, 
                                                        -1, drop = F]
    if (any(is.na(covar))) 
      stop("there are missing values in the covariates")
  }
  ivf <- iv_free(y = y, x = x, g = g, covar = covar, q = q, 
                 family = family)
  x0 <- ivf$x0
  xcoef <- ivf$xcoef
  x0q <- ivf$x0q
  loc <- lace(y = y, x = x, g = g, covar = covar, q = q, x0q = x0q, 
              xc_sub = TRUE, family = family, xpos = xpos)
  coef <- loc$coef/xcoef
  coef_se <- loc$coef_se/xcoef
  xmean <- loc$xmean
  xcoef_sub <- loc$xcoef_sub
  xcoef_sub_se <- loc$xcoef_sub_se
  p_het <- 1 - pchisq(rma(xcoef_sub, vi = (xcoef_sub_se)^2)$QE, 
                      df = (q - 1))
  p_het_trend <- rma.uni(xcoef_sub ~ xmean, vi = xcoef_sub_se^2, 
                         method = method)$pval[2]
  fracpb <- fracpoly_best_powers(coef = coef, coef_se = coef_se, 
                          xmean = xmean, d = d, pd = pd, method = method, powers=powers)
  model <- fracpb$model
  loglik<-fracpb$loglik
  p_ML <- fracpb$p_ML
  p1_ML <- fracpb$p1_ML
  p2_ML <- fracpb$p2_ML
  fp_p <- fracpb$fp_p
  fp_d12_p <- fracpb$fp_d12_p
  d <- fracpb$d
  p_quadratic <- rma(coef ~ xmean, (coef_se)^2, method = "FE")$pval[2]
  p_Q <- 1 - pchisq(rma(coef, vi = (coef_se)^2)$QE, df = (q - 
                                                            1))
  if (ci == "bootstrap_per" | ci == "bootstrap_se") {
    frac_coef_boot <- fracpoly_boot(y = y, x = x, g = g, 
                                    covar = covar, q = q, x0q = x0q, xcoef = xcoef, 
                                    family = family, xpos = xpos, method = method, nboot = nboot, 
                                    d = d, p_ML = p_ML, p1_ML = p1_ML, p2_ML = p2_ML)
  }
  else {
    frac_coef_boot <- NULL
  }
  if (d == 1) {
    powers <- p_ML + 1
  }  else {
    powers <- c((p1_ML + 1), (p2_ML + 1))
  }
  beta <- as.numeric(model$b)
  fitted<-fitted(model)
  if (ci == "model_se") {
    cov <- model$vb
    se <- model$se
    lci <- beta - 1.96 * se
    uci <- beta + 1.96 * se
    pval <- 2 * pnorm(-abs(beta/se))
  }
  if (ci == "bootstrap_se") {
    cov <- var(frac_coef_boot)
    se <- sqrt(diag(cov))
    lci <- beta - 1.96 * se
    uci <- beta + 1.96 * se
    pval <- 2 * pnorm(-abs(beta/se))
  }
  if (ci == "bootstrap_per") {
    if (d == 1) {
      se <- NA
      lci <- quantile(frac_coef_boot, probs = 0.025)
      uci <- quantile(frac_coef_boot, probs = 0.975)
      pval <- NA
    }
    if (d == 2) {
      se <- rep(NA, 2)
      lci <- NULL
      uci <- NULL
      pval <- NULL
      lci[1] <- quantile(frac_coef_boot[, 1], probs = 0.025)
      lci[2] <- quantile(frac_coef_boot[, 2], probs = 0.025)
      uci[1] <- quantile(frac_coef_boot[, 1], probs = 0.975)
      uci[2] <- quantile(frac_coef_boot[, 2], probs = 0.975)
      pval <- rep(NA, 2)
    }
  }
  lci <- as.numeric(lci)
  uci <- as.numeric(uci)
  if (ci == "model_se") {
    nboot <- NA
  }
  if (fig == T) {
    if (d == 1) {
      figure <- fracpoly_figure(beta = beta, cov = cov, 
                                x.min = min(x), x.max = max(x), family = family, 
                                d = d, p_ML = p_ML, ci = ci, frac_coef_boot = frac_coef_boot, 
                                ref = ref, pref_x = pref_x, pref_x_ref = pref_x_ref, 
                                pref_y = pref_y, ci_type = ci_type, ci_quantile = ci_quantile, 
                                breaks = breaks)
    }
    if (d == 2) {
      figure <- fracpoly_figure(beta = beta, cov = cov, 
                                x.min = min(x), x.max = max(x), family = family, 
                                d = d, p1_ML = p1_ML, p2_ML = p2_ML, ci = ci, 
                                frac_coef_boot = frac_coef_boot, ref = ref, 
                                pref_x = pref_x, pref_x_ref = pref_x_ref, pref_y = pref_y, 
                                ci_type = ci_type, ci_quantile = ci_quantile, 
                                breaks = breaks)
    }
  }
  model <- as.matrix(data.frame(q = q, xpos = xpos, ci_type = ci, 
                                nboot = nboot))
  coefficients <- as.matrix(data.frame(beta = beta, se = se, 
                                       lci = lci, uci = uci, pval = pval))
  rownames(coefficients) <- powers
  if (nrow(coefficients) == 2) {
    if (powers[1] == powers[2]) {
      rownames(coefficients) <- c(powers[1], paste0("log ", 
                                                    powers[2]))
    }
  }
  loc <- as.matrix(data.frame(beta = (coef), se = (abs(coef_se)), 
                              lci = (coef - 1.96 * (abs(coef_se))), uci = (coef + 
                                                                             1.96 * (abs(coef_se))), pval = (2 * pnorm(-abs(coef/coef_se)))))
  rownames(loc) <- 1:nrow(loc)
  xcoef_quant <- as.matrix(data.frame(beta = xcoef_sub, se = xcoef_sub_se))
  rownames(xcoef_quant) <- 1:nrow(xcoef_quant)
  p_tests <- as.matrix(data.frame(fp_d1_d2 = fp_d12_p, fp = fp_p, 
                                  quad = p_quadratic, Q = p_Q))
  p_heterogeneity <- as.matrix(data.frame(Q = p_het, trend = p_het_trend))
  if (fig == F) {
    results <- list(n = length(y), model = model, powers = powers, 
                    coefficients = coefficients, lace = loc, xcoef = xcoef_quant, 
                    p_tests = p_tests, p_heterogeneity = p_heterogeneity, loglik=loglik, fitted=fitted)
  }
  else {
    results <- list(n = length(y), model = model, powers = powers, 
                    coefficients = coefficients, lace = loc, xcoef = xcoef_quant, 
                    p_tests = p_tests, p_heterogeneity = p_heterogeneity, loglik=loglik,fitted=fitted,
                    figure = figure)
  }
  class(results) <- "fracpoly_mr"
  return(results)
}
