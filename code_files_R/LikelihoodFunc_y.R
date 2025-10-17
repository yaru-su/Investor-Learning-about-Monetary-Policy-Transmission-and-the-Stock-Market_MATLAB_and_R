# Likelihood function for Output Gap (y)
# Model: dy_t = -λ_y * y_t * Δ + σ_y * sqrt(Δ) * ε_t
# where ε_t ~ N(0,1)

LikelihoodFunc_y <- function(outputgap, param, Delta = 1/12) {
  
  # Number of observations
  n <- length(outputgap)
  
  # Parameters: sigma_y and lambda_y
  sigmay  <- abs(param[1])   # volatility (ensure positive)
  lambday <- param[2]        # mean reversion speed
  
  # Demean output gap (center around zero)
  y <- outputgap - mean(outputgap)
  
  # Initialize storage
  eps <- numeric(n - 1)
  vareps <- numeric(n - 1)
  logLik <- numeric(n - 1)
  stdeps <- numeric(n - 1)
  
  # Loop over time
  for (i in 1:(n - 1)) {
    # Ornstein-Uhlenbeck discrete approximation
    eps[i] <- y[i + 1] - exp(-lambday * Delta) * y[i]
    vareps[i] <- (sigmay^2 / (2 * lambday)) * (1 - exp(-2 * lambday * Delta))
    
    # standardized residuals
    stdeps[i] <- eps[i] / sqrt(vareps[i])
    
    # log-likelihood contribution (Normal)
    Nv <- 1
    logLik[i] <- log(1 / ((2 * pi)^(Nv / 2) * sqrt(vareps[i]))) -
      0.5 * (eps[i]^2 / vareps[i])
  }
  
  # Sum of log-likelihood (negative for minimization)
  sumloglik <- -sum(logLik[-1])  # skip first obs (like MATLAB)
  
  # Return as a list (for flexibility)
  return(list(
    sumloglik = sumloglik,
    logLik = logLik,
    stdeps = stdeps
  ))
}
