# ================================================================
# PURPOSE:
#   Computes Newey-West adjusted heteroskedasticity and serial
#   correlation consistent (HAC) Least Squares regression.
# ================================================================
# USAGE:
#   results <- nwest(y, x, nlag)
# where:
#   y     = dependent variable (nobs x 1 vector)
#   x     = independent variables (nobs x nvar matrix)
#   nlag  = lag length for Newey-West correction
# ================================================================
# RETURNS: a list (same structure as MATLAB version)
#   results$meth           = 'nwest'
#   results$beta           = estimated coefficients
#   results$varianceBeta   = variance of coefficients
#   results$tstat          = t-statistics
#   results$pval           = p-values
#   results$yhat           = fitted values
#   results$resid          = residuals
#   results$sige           = sigma^2
#   results$rsqr           = R^2
#   results$rbar           = adjusted R^2
#   results$dw             = Durbin–Watson statistic
#   results$CovMatBeta     = covariance matrix
#   results$nobs, results$nvar, results$y
# ================================================================

nwest <- function(y, x, nlag = 0) {
  if (missing(y) || missing(x) || missing(nlag)) {
    stop("Usage: results <- nwest(y, x, nlag)")
  }
  
  # Ensure inputs are numeric matrices
  y <- as.matrix(y)
  x <- as.matrix(x)
  nobs <- nrow(x)
  nvar <- ncol(x)
  
  # --------------------------------------------------------------
  # Ordinary Least Squares estimation
  # --------------------------------------------------------------
  xpxi <- solve(t(x) %*% x)
  beta_hat <- xpxi %*% t(x) %*% y
  yhat <- x %*% beta_hat
  resid <- y - yhat
  sigu <- as.numeric(t(resid) %*% resid)
  sige <- sigu / (nobs - nvar)
  
  # --------------------------------------------------------------
  # Construct residual * regressor matrix (hhat)
  # --------------------------------------------------------------
  hhat <- t(resid) %*% x
  hhat <- t(hhat)  # (nvar x nobs)
  
  # --------------------------------------------------------------
  # Newey–West covariance correction
  # --------------------------------------------------------------
  G <- matrix(0, nvar, nvar)
  for (a in 0:nlag) {
    w <- (nlag + 1 - a) / (nlag + 1)   # Bartlett kernel weight
    za <- hhat[, (a + 1):nobs, drop = FALSE] %*%
      t(hhat[, 1:(nobs - a), drop = FALSE])
    if (a == 0) {
      ga <- za
    } else {
      ga <- za + t(za)
    }
    G <- G + w * ga
  }
  
  V <- xpxi %*% G %*% xpxi
  nwerr <- sqrt(diag(V))
  
  # --------------------------------------------------------------
  # Compute t-stats and p-values
  # --------------------------------------------------------------
  tstat <- as.numeric(beta_hat) / nwerr
  pval <- 2 * (1 - pt(abs(tstat), df = nobs - nvar))
  
  # --------------------------------------------------------------
  # R-squared and adjusted R-squared
  # --------------------------------------------------------------
  ym <- y - mean(y)
  rsqr1 <- sigu
  rsqr2 <- as.numeric(t(ym) %*% ym)
  rsqr <- 1 - rsqr1 / rsqr2
  rsqr1 <- rsqr1 / (nobs - nvar)
  rsqr2 <- rsqr2 / (nobs - 1)
  rbar <- 1 - (rsqr1 / rsqr2)
  
  # --------------------------------------------------------------
  # Durbin–Watson statistic
  # --------------------------------------------------------------
  ediff <- resid[2:nobs] - resid[1:(nobs - 1)]
  dw <- as.numeric(t(ediff) %*% ediff / sigu)
  
  # --------------------------------------------------------------
  # Pack results (same structure as MATLAB)
  # --------------------------------------------------------------
  results <- list(
    meth = "nwest",
    beta = as.numeric(beta_hat),
    varianceBeta = nwerr^2,
    tstat = tstat,
    pval = pval,
    yhat = yhat,
    resid = resid,
    sige = sige,
    rsqr = rsqr,
    rbar = rbar,
    dw = dw,
    CovMatBeta = V,
    nobs = nobs,
    nvar = nvar,
    y = y
  )
  
  return(results)
}
