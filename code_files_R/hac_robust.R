# ================================================================
# FUNCTION: hac_robust(reg)
# PURPOSE:
#   Compute Newey–West HAC (Heteroskedasticity and Autocorrelation
#   Consistent) robust standard errors for a linear model and
#   output a formatted table like the MATLAB version.
# ================================================================
# INPUT:
#   reg : lm() regression object
# ================================================================
# OUTPUT:
#   tab_data : data.frame with columns
#     Coeff | SE | tStat | pVal | AdjR2 | Nobs
# ================================================================

hac_robust <- function(reg) {
  # --- 1. Extract residuals ---
  res <- residuals(reg)
  res <- res[!is.na(res)]
  n <- length(res)
  
  # --- 2. Compute autocorrelation and find max significant lag ---
  acf_obj <- acf(res, plot = FALSE)
  acf_vals <- acf_obj$acf[-1]       # skip lag 0
  bounds <- qnorm((1 + 0.95) / 2) / sqrt(n)
  significant <- (acf_vals > bounds) | (acf_vals < -bounds)
  if (any(significant)) {
    max_lag <- which(significant)[length(which(significant))]
  } else {
    max_lag <- 0
  }
  
  # --- 3. Compute HAC robust covariance matrix using QS weights ---
  # Equivalent to MATLAB hac(reg, weights="QS", bandwidth=max_lag)
  cov_hac <- NeweyWest(reg, lag = max_lag, prewhite = FALSE, adjust = TRUE)
  
  # --- 4. Extract coefficients and statistics ---
  coeff <- coef(reg)
  se <- sqrt(diag(cov_hac))
  tstat <- coeff / se
  df <- reg$df.residual
  pval <- 2 * (1 - pt(abs(tstat), df = df))
  
  # --- 5. Construct summary table (same columns as MATLAB) ---
  num_rows <- length(coeff)
  num_cols <- 6
  tab <- matrix(NA, nrow = num_rows, ncol = num_cols)
  
  tab[, 1] <- coeff
  tab[, 2] <- se
  tab[, 3] <- tstat
  tab[, 4] <- pval
  tab[1, 5] <- summary(reg)$adj.r.squared
  tab[1, 6] <- nobs(reg)
  
  colnames(tab) <- c("Coeff", "SE", "tStat", "pVal", "AdjR2", "Nobs")
  tab_data <- as.data.frame(tab)
  
  # --- 6. Return neatly formatted result ---
  return(tab_data)
}
