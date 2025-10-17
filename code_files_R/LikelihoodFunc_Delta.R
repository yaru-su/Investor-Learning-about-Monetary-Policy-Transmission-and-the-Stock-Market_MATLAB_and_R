# Negative log-likelihood for consumption growth
LikelihoodFunc_Delta <- function(ConsGrowth, params){
  mudeltabar <- abs(params[1])
  sigmadelta <- abs(params[2])
  n <- length(ConsGrowth)
  
  sumloglik <- 0
  for(i in 1:(n-1)){
    eps <- ConsGrowth[i] - (mudeltabar - 0.5 * sigmadelta^2)
    sumloglik <- sumloglik - dnorm(eps, mean = 0, sd = sigmadelta, log = TRUE)
  }
  
  return(sumloglik)
}