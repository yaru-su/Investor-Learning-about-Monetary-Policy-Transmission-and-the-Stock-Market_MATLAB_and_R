GetVarMatrixParam_Rate <- function(fedfund, inflation, outputgap, param) {
  # --------------------------------------------------------------
  # GetVarMatrixParam_Rate:
  # Computes the variance-covariance matrix of the estimated parameters
  # using the outer product of the gradient (score) method.
  #
  # Methods implemented (only method 2 used):
  # (1) Second partial derivatives (not implemented here)
  # (2) First derivatives (outer product matrix) [used]
  # (3) White covariance matrix (not used)
  # (4) Newey-West covariance matrix (not used)
  # --------------------------------------------------------------
  
  # Small perturbation for finite difference derivative approximation
  myDelta <- 1e-15 * abs(param)
  
  # Baseline likelihood and log-likelihood vector
  base <- LikelihoodFunc_Rate(fedfund, inflation, outputgap, param)
  logLikVec1 <- base$logLik
  n <- length(logLikVec1)
  
  # Initialize score matrix (each column = observation, each row = parameter derivative)
  s <- matrix(0, nrow = length(param), ncol = n)
  
  # Finite difference: compute numerical derivative for each parameter
  for (i in seq_along(param)) {
    vec <- param
    vec[i] <- param[i] + myDelta[i]
    
    out <- LikelihoodFunc_Rate(fedfund, inflation, outputgap, vec)
    logLikVec2 <- out$logLik
    
    s[i, ] <- (logLikVec2 - logLikVec1) / myDelta[i]
  }
  
  # Outer product of gradients
  OP_Matrix <- (s %*% t(s)) / n
  
  # Inverse of OPG (variance-covariance matrix)
  V <- solve(OP_Matrix) / n
  
  return(V)
}
