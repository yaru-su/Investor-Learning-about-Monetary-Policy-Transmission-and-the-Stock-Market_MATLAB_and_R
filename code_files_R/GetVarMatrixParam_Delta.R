# Variance-covariance matrix of estimated parameters
GetVarMatrixParam_Delta <- function(ConsGrowth, params){
  # Compute Hessian
  hess <- hessian(func = LikelihoodFunc_Delta, x = params, ConsGrowth = ConsGrowth)
  # Invert Hessian to get covariance matrix
  V <- solve(hess)
  return(V)
}