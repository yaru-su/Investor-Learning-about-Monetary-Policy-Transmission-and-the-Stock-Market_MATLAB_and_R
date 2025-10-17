LikelihoodFunc_Rate <- function(fedfund, inflation, outputgap, param) {
  # ---------------------------------------------------------------
  # LikelihoodFunc_Rate: Log-likelihood for Taylor Rule estimation
  # Input:
  #   fedfund   - observed nominal interest rate
  #   inflation - inflation rate series
  #   outputgap - output gap series
  #   param     - parameter vector (rNbar, betapi, betay, volr)
  # Output:
  #   sumloglik - total negative log-likelihood (for optimization)
  #   logLik    - vector of log-likelihood values
  #   rfTaylor  - Taylor rule implied interest rate
  #   stdeps    - standardized residuals
  # ---------------------------------------------------------------
  
  n <- length(fedfund)
  
  # Extract parameters
  rNbar <- param[1]
  betapi <- abs(param[2])
  betay <- abs(param[3])
  volr <- abs(param[4])   # standard deviation of residuals (not used directly in model)
  
  meaninflation <- mean(inflation)
  meanoutputgap <- mean(outputgap)
  
  # Initialize output vectors
  rfTaylor <- numeric(n)
  eps <- numeric(n)
  stdeps <- numeric(n)
  logLik <- numeric(n)
  
  for (i in seq_len(n)) {
    # Compute model-implied interest rate
    rfTaylor[i] <- rTaylor(
      rNbar, betapi, betay,
      inflation[i], outputgap[i],
      meaninflation, meanoutputgap
    )
    
    # Residuals
    eps[i] <- fedfund[i] - rfTaylor[i]
    vareps <- volr^2
    stdeps[i] <- eps[i] / sqrt(vareps)
    
    # Normal log-likelihood for each observation
    Nv <- 1
    logLik[i] <- log(1 / ((2 * pi)^(Nv / 2) * sqrt(vareps))) - 0.5 * (eps[i]^2 / vareps)
  }
  
  # Apply penalty if model misbehaves (same as MATLAB version)
  if (abs(mean(rfTaylor) - mean(fedfund)) > 0.00005 ||
      sd(betapi * inflation) > 4 * sd(betay * outputgap)) {
    pen <- 1e10
  } else {
    pen <- 0
  }
  
  # Sum of log-likelihood (negative for minimization)
  sumloglik <- -sum(logLik[-1]) + pen
  
  # Return results as list
  return(list(
    sumloglik = sumloglik,
    logLik = logLik,
    rfTaylor = rfTaylor,
    stdeps = stdeps
  ))
}
