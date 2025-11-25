# -------------------------------------------------------------
# Function: GetVarMatrixParamy
# Purpose : Compute the variance-covariance matrix of parameters 
#           using the outer-product of likelihood gradients.
# -------------------------------------------------------------
GetVarMatrixParam_y <- function(outputgap, param) {
  # Small perturbation for finite differences (two-sided derivative)
  myDelta <- 1e-15 * abs(param)
  
  # Evaluate log-likelihood at the current parameter
  ll_current <- LikelihoodFunc_y(outputgap, param)
  logLikVec1 <- ll_current$logLik  # vector of log-likelihoods per observation
  
  n <- length(logLikVec1)
  k <- length(param)
  s <- matrix(0, nrow = k, ncol = n)  # store score vectors
  
  # Loop over parameters for finite-difference derivative
  for (i in seq_len(k)) {
    vec <- param
    vec[i] <- param[i] + myDelta[i]  # perturb parameter i
    
    ll_perturbed <- LikelihoodFunc_y(outputgap, vec)
    logLikVec2 <- ll_perturbed$logLik
    
    # Compute numerical derivative (first derivative of log-likelihood)
    s[i, ] <- (logLikVec2 - logLikVec1) / myDelta[i]
  }
  
  # Outer-product matrix of scores (average information matrix)
  OP_Matrix <- (s %*% t(s)) / n
  
  # Variance-covariance matrix of parameters
  V <- solve(OP_Matrix) / n
  
  return(V)
}
