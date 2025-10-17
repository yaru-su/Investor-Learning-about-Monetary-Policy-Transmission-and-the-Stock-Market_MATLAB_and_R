GetVarMatrixParam <- function(inflation, fedfund, param) {
  # --------------------------------------------------------------
  # GetVarMatrixParam:
  #   Computes the variance-covariance matrix of estimated parameters
  #   using a two-sided finite difference approximation of derivatives.
  #
  # Four covariance estimation methods (only method 2 used here):
  #   (1) Second partial derivatives
  #   (2) First partial derivatives (Outer Product Matrix)
  #   (3) White covariance matrix (NOT USED)
  #   (4) Newey-West covariance matrix (NOT USED)
  #
  # References:
  #   Hamilton, J.D. (1994). "Time Series Analysis", Princeton University Press.
  #   Newey, W. & West, K. (1987). "A Simple, Positive Semidefinite,
  #   Heteroskedasticity and Autocorrelation Consistent Covariance Matrix",
  #   Econometrica, 55, 347–370.
  #   White, H. (1984). "Asymptotic Theory for Econometricians", Academic Press.
  # --------------------------------------------------------------
  
  # Small perturbation for finite difference derivative approximation
  myDelta <- 1e-15 * abs(param)
  
  # First derivative calculation
  base <- LikelihoodFunc(inflation, fedfund, param)
  logLikVec1 <- base$logLik
  n <- length(logLikVec1)
  s <- matrix(0, nrow = length(param), ncol = n)
  
  # Compute numerical score (gradient) for each parameter by finite difference
  for (i in seq_along(param)) {
    vec <- param
    vec[i] <- param[i] + myDelta[i]
    
    out <- LikelihoodFunc(inflation, fedfund, vec)
    logLikVec2 <- out$logLik
    
    s[i, ] <- (logLikVec2 - logLikVec1) / myDelta[i]
  }
  
  # Outer product of gradients (score matrix)
  OP_Matrix <- (s %*% t(s)) / n
  
  # Invert the outer product matrix to obtain covariance matrix
  V <- solve(OP_Matrix) / n
  
  return(V)
}
