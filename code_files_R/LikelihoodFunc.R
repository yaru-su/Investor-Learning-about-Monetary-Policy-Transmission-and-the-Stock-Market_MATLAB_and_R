LikelihoodFunc <- function(inflation, fedfunds, param) {
  # --------------------------------------------------------------
  # LikelihoodFunc:
  #   Computes the (negative) log-likelihood for the inflation block
  #   based on an Ornstein-Uhlenbeck process and a Kalman learning model.
  #
  # Inputs:
  #   inflation - observed inflation time series
  #   fedfunds  - Federal funds rate series
  #   param     - parameter vector (sigmapi, pibar, lambdapi, sigmaa, lambdaa)
  #
  # Output (list):
  #   sumloglik - negative total log-likelihood (scalar)
  #   logLik    - vector of log-likelihood values per observation
  #   a         - learning coefficient sequence
  #   nua       - uncertainty sequence
  #   LTMpifinal - long-term mean of inflation
  #   stdeps    - standardized residuals
  # --------------------------------------------------------------
  
  # Use global Delta (time step)
  Delta <- get("Delta", envir = .GlobalEnv)
  
  n <- length(inflation)
  
  # Extract parameters (ensure nonnegativity where required)
  sigmapi <- abs(param[1])
  pibar   <- param[2]
  lambdapi <- abs(param[3])
  sigmaa  <- abs(param[4])
  lambdaa <- abs(param[5])
  abar <- 0
  rNbar <- mean(fedfunds)
  
  # Initialize state variables
  a <- numeric(n)
  nua <- numeric(n)
  a[1] <- abar
  nua[1] <- 0.5 * sigmaa^2 / (2 * lambdaa)
  
  # Initialize storage vectors
  LTMpi <- numeric(n - 1)
  eps <- numeric(n - 1)
  stdeps <- numeric(n - 1)
  logLik <- numeric(n - 1)
  
  # Main likelihood loop
  for (i in 1:(n - 1)) {
    # Approximate solution of Ornstein-Uhlenbeck (OU)
    LTMpi[i] <- pibar - a[i] * (fedfunds[i] - rNbar)
    eps[i] <- inflation[i + 1] -
      exp(-lambdapi * Delta) * inflation[i] -
      LTMpi[i] * (1 - exp(-lambdapi * Delta))
    vareps <- sigmapi^2 / (2 * lambdapi) * (1 - exp(-2 * lambdapi * Delta))
    stdeps[i] <- eps[i] / sqrt(vareps)
    
    # Update a and nua
    vara <- ((rNbar - fedfunds[i]) * lambdapi * nua[i] / sigmapi)^2 /
      (2 * lambdaa) * (1 - exp(-2 * lambdaa * Delta))
    a[i + 1] <- exp(-lambdaa * Delta) * a[i] +
      abar * (1 - exp(-lambdaa * Delta)) +
      sign(rNbar - fedfunds[i]) * sqrt(vara) * stdeps[i]
    
    nua[i + 1] <- max(
      0,
      nua[i] + (sigmaa^2 -
                  2 * lambdaa * nua[i] -
                  (rNbar - fedfunds[i])^2 * lambdapi^2 * nua[i]^2 / sigmapi^2) * Delta
    )
    
    # Normal log-likelihood for this observation
    logLik[i] <- log(1 / (sqrt(2 * pi * vareps))) - 0.5 * (eps[i]^2 / vareps)
  }
  
  # Penalty conditions (same logic as MATLAB)
  if (
    mean(LTMpi) + 2 * sd(LTMpi) > quantile(inflation, 0.9) ||
    mean(LTMpi) - 2 * sd(LTMpi) < quantile(inflation, 0.1) ||
    mean(LTMpi) > mean(inflation) + 0.0005 ||
    mean(LTMpi) < mean(inflation) - 0.0005 ||
    lambdaa < 0.5
  ) {
    pen <- 1e10
  } else {
    pen <- 0
  }
  
  # Sum of log-likelihood (negative for optimization)
  sumloglik <- -sum(logLik[-1]) + pen
  
  # Compute final long-term inflation mean
  LTMpifinal <- pibar - a * (fedfunds - rNbar)
  
  # Return list
  return(list(
    sumloglik = sumloglik,
    logLik = logLik,
    a = a,
    nua = nua,
    LTMpifinal = LTMpifinal,
    stdeps = stdeps
  ))
}
