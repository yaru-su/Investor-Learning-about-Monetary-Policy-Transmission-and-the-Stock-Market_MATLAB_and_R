library(tidyverse)
library(lubridate)
library(patchwork)
library(numDeriv)
library(stats)  # for optim, pnorm, qt functions
library(sandwich)   # For Newey-West covariance
library(lmtest)     # For coeftest()

setwd("/Users/yarusu/Library/CloudStorage/OneDrive-UCSanDiego/DS/ML project/AndreiHaslerJFEReplication_V2/MatlabCodeAndData")

# ==============================================================================
# 1. Figure 1
# ==============================================================================

# Convert data to long format for faceted plotting
TT_long <- TT %>%
  select(date, GDPGrowth, OutputGap, inflation, rfnom) %>%
  pivot_longer(
    cols = -date,                 # Keep 'date' as the identifier
    names_to = "variable",        # Column names become a new 'variable' column
    values_to = "value"           # Values go into a new 'value' column
  )

# Assign y-axis labels for each variable
y_labels <- c(
  GDPGrowth = "GDP growth",
  OutputGap = "Output gap",
  inflation = "Inflation",
  rfnom = "Fed funds rate"
)

# Create the faceted line plot
p <- ggplot(TT_long, aes(x = date, y = value)) +
  geom_line(color = "blue") +                 # Draw line in blue
  facet_wrap(~variable, ncol = 1, scales = "free_y", 
             labeller = labeller(variable = y_labels)) +  # Label panels
  labs(title = "Figure 1: GDP Growth, Output Gap, Inflation, and Federal Funds Rate") +
  theme_minimal(base_size = 14)              # Clean minimal theme

# Save the figure as a PNG
ggsave("Figure1_MacroVariables_R.png", p, width = 10, height = 8)


# ==============================================================================
# 2. Consumption Growth: ML Estimation of (mudeltabar, sigmadelta) mudeltabar not needed in the model
# ==============================================================================
# Requires: LikelihoodFunc_Delta_R, GetVarMatrixParam_Delta

# Initial guess for parameters [mudeltabar, sigmadelta]
paramdelta0 <- c(0.02, 0.02)

# Maximize likelihood
fit <- optim(
  par = paramdelta0,
  fn = LikelihoodFunc_Delta,
  ConsGrowth = TT$GDPGrowth,
  method = "L-BFGS-B",
  control = list(fnscale = 1, reltol = 1e-4)
)

# Estimated parameters
paramdelta <- fit$par

# Compute variance-covariance matrix
Vdelta <- GetVarMatrixParam_Delta(TT$GDPGrowth, paramdelta)

# Standard errors
paramdelta_std <- sqrt(diag(Vdelta))

# Number of observations
nobs <- length(TT$GDPGrowth)

# T-statistics
paramdelta_tstats <- abs(paramdelta / paramdelta_std)

# Two-sided p-values
paramdelta_pvalues <- 2 * (1 - pt(paramdelta_tstats, df = nobs - length(paramdelta)))

# Display results
cat("Estimated Parameters:\n"); print(paramdelta)
cat("Standard Errors:\n"); print(paramdelta_std)
cat("T-Statistics:\n"); print(paramdelta_tstats)
cat("P-Values:\n"); print(paramdelta_pvalues)


# ==============================================================================
# 3. Output Gap: ML Estimation of (mudeltabar, sigmadelta) mudeltabar not needed in the model
# ==============================================================================
# Requires: LikelihoodFunc_y, GetVarMatrixParam_y

# Estimate initial guess for lambda_y using AR(1) regression
y <- TT$OutputGap
y_lag <- head(y, -1)
y_lead <- tail(y, -1)

# Simple OLS regression of y_t on y_{t-1}
reg_model <- lm(y_lead ~ y_lag)
reg_coef <- coef(reg_model)[2]   # AR(1) coefficient

# Convert AR(1) coefficient to continuous-time mean reversion rate
lamby_reg <- -log(reg_coef) / Delta

# Set initial parameter vector (sigma_y, lambda_y)
paramy0 <- c(sd(y, na.rm = TRUE), lamby_reg)

# Maximize likelihood function to estimate parameters
obj_fun <- function(par) {
  LikelihoodFunc_y(y, par)$sumloglik   # Returns negative log-likelihood
}

# Use Nelder-Mead
optim_result <- optim(paramy0, obj_fun, method = "Nelder-Mead")
paramy <- optim_result$par
fval <- optim_result$value

# Compute variance-covariance matrix of estimated parameters
Vy <- GetVarMatrixParam_y(y, paramy)

# Compute standard errors, t-stats, and p-values
paramy_std <- sqrt(diag(Vy))
nobs <- length(y)
paramy_tstats <- abs(paramy / paramy_std)
paramy_pvalues <- 2 * (1 - pt(paramy_tstats, df = nobs - length(paramy)))

# Display results
cat("Estimated Parameters for Output Gap:\n")
print(paramy)
cat("Standard Errors for Output Gap:\n")
print(paramy_std)
cat("T-Statistics for Output Gap:\n")
print(paramy_tstats)
cat("P-Values for Output Gap:\n")
print(paramy_pvalues)


# ==============================================================================
# 4. Taylor Rule: ML Estimation of (rNbar,betapi,betay,volr) volr not needed in model
# ==============================================================================
# Requires: LikelihoodFunc_Rate, GetVarMatrixParam_Rate

# Define initial parameter guesses
paramr0 <- c(mean(TT$rfnom), 0.25, 0.15, 0.02)

# Define the negative log-likelihood function for optimization
negLogLik <- function(x) {
  LikelihoodFunc_Rate(TT$rfnom, TT$inflation, TT$OutputGap, x)
}

# Use optim (similar to fminsearch) to perform Maximum Likelihood estimation
result <- optim(
  par = paramr0,
  fn = function(x) {
    out <- LikelihoodFunc_Rate(TT$rfnom, TT$inflation, TT$OutputGap, x)
    return(out$sumloglik)
  },
  method = "Nelder-Mead"
)

# Extract estimated parameters
paramr <- result$par
fval <- result$value

# Compute log-likelihood and model-implied variables
out <- LikelihoodFunc_Rate(TT$rfnom, TT$inflation, TT$OutputGap, paramr)
sumloglik <- out$sumloglik
logLik <- out$logLik
rfTaylor <- out$rfTaylor
stdeps <- out$stdeps

# Compute variance-covariance matrix of estimated parameters
Vr <- GetVarMatrixParam_Rate(TT$rfnom, TT$inflation, TT$OutputGap, paramr)

# Compute standard errors from the diagonal of the variance-covariance matrix
paramr_std <- sqrt(diag(Vr))

# Number of observations
nobs <- length(TT$inflation)

# Compute t-statistics for estimated parameters
paramr_tstats <- abs(paramr / paramr_std)

# Compute two-sided p-values
paramr_pvalues <- 2 * (1 - pt(paramr_tstats, df = nobs - length(paramr)))

# Extract estimated rNbar
rNbar_estimated <- paramr[1]

# Display results
cat("Estimated Parameters for Taylor Rule:\n")
print(paramr)
cat("Standard Errors for Taylor Rule:\n")
print(paramr_std)
cat("T-Statistics for Taylor Rule:\n")
print(paramr_tstats)
cat("P-Values for Taylor Rule:\n")
print(paramr_pvalues)


# ==============================================================================
# 5. Inflation and Transmission Mechanism: ML estimation of inflation parameters (sigmapi,pibar,lambdapi,sigmaa,lambdaa)
# ==============================================================================
# Requires: LikelihoodFunc, GetVarMatrixParam

# Compute sample covariance between current and lagged inflation
COV <- cov(cbind(TT$inflation[-1], TT$inflation[-length(TT$inflation)]))

# Compute initial mean-reversion speed λπ
lambdapi0 <- -1 / Delta * log(COV[1, 2] / COV[1, 1])

# Compute initial volatility σπ
sigmapi0 <- sqrt(var(TT$inflation) * 2 * lambdapi0)

# Compute initial long-run inflation target π̄
pibar0 <- mean(TT$inflation)

# Combine all initial guesses into parameter vector
param0 <- c(sigmapi0, pibar0, 0.1, 0.2, 0.5)

# Estimate parameters via Maximum Likelihood
result <- optim(
  par = param0,
  fn = function(x) {
    out <- LikelihoodFunc(TT$inflation, TT$fedfund, x)
    return(out$sumloglik)
  },
  method = "Nelder-Mead"
)

param <- result$par
fval <- result$value

# Recompute likelihood and filtered states using estimated parameters
out <- LikelihoodFunc(TT$inflation, TT$fedfund, param)
sumloglik <- out$sumloglik
logLik <- out$logLik
a <- out$a
nua <- out$nua
LTMpi <- out$LTMpi
stdeps <- out$stdeps

# Compute variance-covariance matrix of parameters
V <- GetVarMatrixParam(TT$inflation, TT$fedfund, param)

# Compute standard errors from the diagonal of the variance-covariance matrix
param_std <- sqrt(diag(V))

# Number of observations
nobs <- length(TT$inflation)

# Compute t-stat of estimated parameters
param_tstats <- abs(param / param_std)

# Compute two-sided p-values
param_pvalues <- 2 * (1 - pt(param_tstats, df = nobs - length(param)))

# Display estimated parameters and their statistics for inflation
cat("Estimated Parameters for Inflation:\n")
print(param)
cat("Standard Errors for Inflation:\n")
print(param_std)
cat("T-Statistics for Inflation:\n")
print(param_tstats)
cat("P-Values for Inflation:\n")
print(param_pvalues)

# Insert ML-Kalman estimated state variables in timetable
TT$a <- a                # learning (transmission) coefficient
TT$nua <- nua            # uncertainty about a_t
TT$LTMpi <- LTMpi        # long-term mean of inflation
TT$phi <- rfTaylor - rNbar_estimated   # monetary policy deviation

# Analyze properties of estimated state variables
# (a) Relationship between changes in a_t and inflation
change_a <- diff(TT$a)
change_pi <- diff(TT$inflation)
dummy <- as.numeric(TT$inflation > mean(TT$inflation))
corr_high <- cor(change_a[dummy[-1] == 1], change_pi[dummy[-1] == 1])
corr_low <- cor(change_a[dummy[-1] == 0], change_pi[dummy[-1] == 0])

# (b) Moments of φ = r_N,t - r̄_N (monetary policy deviation)
phi_lead <- TT$phi[-1]
phi_lag <- TT$phi[-length(TT$phi)]
reg_autocorr_phi <- lm(phi_lead ~ phi_lag)
autocorrphi <- coef(reg_autocorr_phi)[2]
Volphi <- sd(TT$phi, na.rm = TRUE)

# (c) Moments of long-term mean of inflation (LTMπ)
LTMpi_lead <- LTMpi[-1]
LTMpi_lag <- LTMpi[-length(LTMpi)]
reg_LTMpi <- lm(LTMpi_lead ~ LTMpi_lag)
autocorLTMpi <- coef(reg_LTMpi)[2]
volLTMpi <- sd(LTMpi, na.rm = TRUE)

# (d) Correlation with 5-year inflation expectations (footnote 11)
corr_CPI5YR <- cor(TT$CPI5YR, TT$LTMpi, use = "complete.obs")
reg_CPI5YR <- lm(CPI5YR ~ LTMpi, data = TT)

# HAC robust regression (Newey-West)
library(sandwich)
library(lmtest)
hac_results <- coeftest(reg_CPI5YR, vcov = NeweyWest(reg_CPI5YR, prewhite = FALSE))

# Display analysis results
cat("Correlation (high inflation periods):", corr_high, "\n")
cat("Correlation (low inflation periods):", corr_low, "\n")
cat("Autocorrelation of phi:", autocorrphi, "\n")
cat("Std. deviation of phi:", Volphi, "\n")
cat("Autocorrelation of LTMpi:", autocorLTMpi, "\n")
cat("Std. deviation of LTMpi:", volLTMpi, "\n")
cat("Correlation between CPI5YR and LTMpi:", corr_CPI5YR, "\n")
cat("HAC-robust regression results (CPI5YR ~ LTMpi):\n")
print(hac_results)


# ==============================================================================
# 6. Figure 2: Investors' Estimate of the Transmission Coefficient
# ==============================================================================

# Save the figure as a PNG
png("Figure2_TransmissionCoefficient_R.png", width = 10, height = 8, units = "in", res = 300)

# Left y-axis: transmission coefficient
plot(TT$date, TT$a,
     type = "l", col = "blue", lwd = 2,
     ylab = expression("Transmission coefficient " * hat(a)[t]),
     xlab = "Date")

# Right y-axis: inflation
par(new = TRUE)  # overlay new plot on existing figure
plot(TT$date, TT$inflation,
     type = "l", col = "red", lwd = 2, lty = 2,
     axes = FALSE, xlab = "", ylab = "")

axis(side = 4, col.axis = "red", col = "red")
mtext(expression("Inflation " * pi[t]), side = 4, line = 3, col = "red")

# Add legend
legend("bottomleft",
       legend = c(expression("Transmission coefficient " * hat(a)[t]),
                  expression("Inflation " * pi[t])),
       col = c("blue", "red"),
       lwd = 2, lty = c(1, 2),
       bty = "n")

# Close PNG device (save file)
dev.off()


# ==============================================================================
# 7. Figure 3: Monetary Policy Stance and Long-Term Inflation Drift
# ==============================================================================

# Save the figure as a PNG
png("Figure3_MonetaryPolicy_LTMInflation_R.png", width = 10, height = 8, units = "in", res = 300)

# Set margins and layout (2 rows, 1 column)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2) + 0.1)

# First panel
plot(TT$date, TT$phi,
     type = "l", col = "blue", lwd = 2,
     ylab = expression("Monetary policy stance " * phi[t]),
     xlab = "Date")
abline(h = 0, col = "gray60", lty = 2)  # reference line at φ_t = 0

# Second panel
plot(TT$date, TT$LTMpi,
     type = "l", col = "red", lwd = 2,
     ylab = expression("Long-term inflation drift " * hat(pi)[t]^A),
     xlab = "Date")
abline(h = mean(TT$LTMpi, na.rm = TRUE), col = "gray60", lty = 2)  # mean reference line

# Close PNG device (save file)
dev.off()


# ==============================================================================
# 8. Table 3: Parameter Values Estimated by Maximum Likelihood
# ==============================================================================

# Columns: sigmadelta, sigmay, lambday, rNbar, betapi, betay, sigmapi, pibar, lambdapi, sigmaa, lambdaa

# Combine estimated parameters and standard errors into a summary table
ParamEstimTable <- rbind(
  c(paramdelta[2], paramdelta_std[2]),
  cbind(paramy, paramy_std),
  cbind(paramr[1:(length(paramr) - 1)], paramr_std[1:(length(paramr_std) - 1)]),
  cbind(param, param_std)
)

# Assign column and row names for readability
colnames(ParamEstimTable) <- c("Estimate", "Std.Error")
rownames(ParamEstimTable) <- c(
  "sigmadelta",
  "sigmay", "lambday",
  "rNbar", "betapi", "betay",
  "sigmapi", "pibar", "lambdapi", "sigmaa", "lambdaa"
)

# Display parameter summary table (equivalent to MATLAB variable ParamEstimTable)
cat("Parameter Estimation Table (Table 3):\n")
print(round(ParamEstimTable, 5))


# ==============================================================================
# 9. SAVE ML-KALMAN ESTIMATED STATE VARIABLES
# ==============================================================================
StateVars <- TT[, c("GDPGrowth", "rfnom", "inflation", "a", "nua", "phi", "OutputGap", "LTMpi")]

# Export to CSV (same structure as MATLAB writetimetable)
write.csv(StateVars, "StateVars_Oct2024.csv", row.names = FALSE)
cat("\nState variables successfully saved to 'StateVars_Oct2024.csv'\n")


# ==============================================================================
# 10. S&P500 real yearly log dividend growth (used to calibrate sigmaD in model)
# ==============================================================================

# Compute real (inflation-adjusted) yearly log dividend growth
divgrowth <- log(TT$d12[13:length(TT$d12)] / TT$d12[1:(length(TT$d12) - 12)]) -
  TT$inflation[13:length(TT$inflation)]

# Compute standard deviation of dividend growth (σ_D)
sigmaD <- sd(divgrowth, na.rm = TRUE)

# Display result
cat("\nStandard deviation of real yearly log dividend growth (sigmaD):", round(sigmaD, 6), "\n")


# ==============================================================================
# 11. Table 4: Asset-pricing Moments
# ==============================================================================

# (a) Data Summary Statistics

# Real and nominal risk-free rates
meanrfreal <- mean(TT$rfreal, na.rm = TRUE)
meanrfnom  <- mean(TT$rfnom,  na.rm = TRUE)

# S&P500 excess returns (annualized)
meanSP <- 12 * mean(TT$exret, na.rm = TRUE)
volSP  <- sqrt(12) * sd(TT$exret, na.rm = TRUE)

# Display results
cat("Table 4: Data Moments\n")
cat("Mean real risk-free rate:", meanrfreal, "\n")
cat("Mean nominal risk-free rate:", meanrfnom, "\n")
cat("Mean annualized S&P500 excess return:", meanSP, "\n")
cat("Annualized volatility of S&P500 excess return:", volSP, "\n")

# (b) Model-Implied Time Series (from Mathematica output)

# Load model-implied data (from Mathematica export)
DataModel <- read.csv("data/ModelImpliedTimeSeriesNov2024.csv")

# Align model-implied data with empirical sample dates
DataModel$date <- TT$date   # ensure same timeline as empirical TT

# Convert to time series tibble (optional)
library(tibbletime)
TTmodel <- as_tibble(DataModel)

# Display confirmation
cat("Model-implied time series loaded and aligned with empirical data.\n")
cat("Number of observations:", nrow(TTmodel), "\n")
cat("Variables:", paste(colnames(TTmodel), collapse = ", "), "\n")

# ==============================================================================
# 12. Table 5: Expected Output Growth, Real Interest Rate, and Log Price-Dividend Ratio (X) vs. Inflation and the Output Gap (Y).
# ==============================================================================

## -----------------------------------------------------------------------------
## (a) Col 1,2: Expected Output Growth
## -----------------------------------------------------------------------------

### Model regression
# Run OLS regression: regress model-implied expected output growth on inflation and output gap.
reg <- lm(TTmodel$mudelta ~ TT$inflation + TT$OutputGap)

# Compute HAC-robust (Newey–West) standard errors for time series data.
tab_model <- hac_robust(reg)

# Identify which observations (rows) were used (exclude missing values).
usedObs <- complete.cases(TT$inflation, TT$OutputGap, TTmodel$mudelta)

# Extract the corresponding dates from TT
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Compute standardized effect of a one-standard-deviation increase in inflation and output gap on expected output growth.
OneStdIncrease_x <- apply(cbind(TT$inflation[usedObs], TT$OutputGap[usedObs]), 2, sd, na.rm = TRUE) *
  coef(reg)[2:3]
OneStd_y <- sd(TTmodel$mudelta[usedObs], na.rm = TRUE)
OneStdIncrease_x / OneStd_y


###  Data regression
reg <- lm(TT$mudelta ~ TT$inflation + TT$OutputGap)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$inflation, TT$OutputGap, TT$mudelta)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)
OneStdIncrease_x <- apply(cbind(TT$inflation[usedObs], TT$OutputGap[usedObs]), 2, sd, na.rm = TRUE) *
  coef(reg)[2:3]
OneStd_y <- sd(TT$mudelta[usedObs], na.rm = TRUE)
OneStdIncrease_x / OneStd_y

# Table 5, col 1,2
tab_model
tab_data


## -----------------------------------------------------------------------------
## (b) Col 3,4: Real Risk-Free Rate
## -----------------------------------------------------------------------------

### Model regression
reg <- lm(TTmodel$rfreal ~ TT$inflation + TT$OutputGap)
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$inflation, TT$OutputGap, TTmodel$rfreal)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)
OneStdIncrease_x <- apply(cbind(TT$inflation[usedObs], TT$OutputGap[usedObs]), 2, sd, na.rm = TRUE) *
  coef(reg)[2:3]
OneStd_y <- sd(TTmodel$rfreal[usedObs], na.rm = TRUE)
OneStdIncrease_x / OneStd_y


### Data regression
reg <- lm(TT$rfreal ~ TT$inflation + TT$OutputGap)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$inflation, TT$OutputGap, TT$rfreal)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)
OneStdIncrease_x <- apply(cbind(TT$inflation[usedObs], TT$OutputGap[usedObs]), 2, sd, na.rm = TRUE) *
  coef(reg)[2:3]
OneStd_y <- sd(TT$rfreal[usedObs], na.rm = TRUE)
OneStdIncrease_x / OneStd_y

# Table 5, col 3,4
tab_model
tab_data


## -----------------------------------------------------------------------------
## (c) Col 5,6: Log PD Ratio
## -----------------------------------------------------------------------------

### Model regression
reg <- lm(TTmodel$pd ~ TT$inflation + TT$OutputGap)
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$inflation, TT$OutputGap, TTmodel$pd)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)
OneStdIncrease_x <- apply(cbind(TT$inflation[usedObs], TT$OutputGap[usedObs]), 2, sd, na.rm = TRUE) *
  coef(reg)[2:3]
OneStd_y <- sd(TTmodel$pd[usedObs], na.rm = TRUE)
OneStdIncrease_x / OneStd_y


### Data regression
reg <- lm(TT$pd ~ TT$inflation + TT$OutputGap)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$inflation, TT$OutputGap, TT$pd)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)
OneStdIncrease_x <- apply(cbind(TT$inflation[usedObs], TT$OutputGap[usedObs]), 2, sd, na.rm = TRUE) *
  coef(reg)[2:3]
OneStd_y <- sd(TT$pd[usedObs], na.rm = TRUE)
OneStdIncrease_x / OneStd_y

# Table 5, col 5,6
tab_model
tab_data


# ==============================================================================
# 13. Table 6:  Market Risk Premium and Return Volatility (X) vs. State Variables (Y)
# ==============================================================================

## -----------------------------------------------------------------------------
## Panel A
## -----------------------------------------------------------------------------

### Panel A, col 1,2: Model regression
### ------------------------------------
# corr between Martin(model) and ChabiYo RPs(data)
# 'use = "complete.obs"' means ignore rows with missing values (NaN).
cor(TT$RP1, 12 * TT$LBR_30, use = "complete.obs")

# Define risk premium variables
RPmodel <- TTmodel$RP           # Model-implied risk premium from theoretical model (TTmodel)
RPemp   <- 12 * TT$LBR_30       # Annualized risk premium (1-month-ahead, preference-free) ChabiYo

# Fit a linear model (OLS): RPmodel_t = α + β * a_t + ε_t
reg <- lm(RPmodel ~ TT$a)

# Compute HAC-robust (Newey–West) standard errors
tab_model <- hac_robust(reg)

# Identify which observations (rows) were used (exclude missing values)
usedObs <- complete.cases(TT$a, RPmodel)

# Extract the corresponding dates from TT
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)


### Panel A, col 1,2: Data regression
### ------------------------------------
reg <- lm(RPemp ~ TT$a)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$a, RPemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel A, col 1,2
tab_model
tab_data


### Panel A, col 3,4: Model regression
### ------------------------------------
reg <- lm(RPmodel ~ I(TT$phi^2))
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$phi, RPmodel)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

### Panel A, col 3,4: Data regression
### ------------------------------------
reg <- lm(RPemp ~ I(TT$phi^2))
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$phi, RPemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel A, col 3,4
tab_model
tab_data


### Panel A, col 5,6: Model regression
### ------------------------------------
reg <- lm(RPmodel ~ TT$a + I(TT$phi^2))
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, RPmodel)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

### Panel A, col 5,6: Data regression
### ------------------------------------
reg <- lm(RPemp ~ TT$a + I(TT$phi^2))
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, RPemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel A, col 5,6
tab_model
tab_data


### Panel A, col 7,8: Model regression
### ------------------------------------
reg <- lm(RPmodel ~ TT$a + TT$phi + I(TT$phi^2))
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, RPmodel)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Compute standardized effect of a one-standard-deviation increase in a, phi, and phi^2
Xmat <- cbind(TT$a[usedObs], TT$phi[usedObs], TT$phi[usedObs]^2)
OneStdIncrease_x <- apply(Xmat, 2, sd, na.rm = TRUE) * coef(reg)[2:4]
OneStd_y <- sd(RPmodel[usedObs], na.rm = TRUE)
OneStdIncrease_x / OneStd_y

### Panel A, col 7,8: Data regression
### ------------------------------------
reg <- lm(RPemp ~ TT$a + TT$phi + I(TT$phi^2))
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, RPemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Compute standardized effect of a one-standard-deviation increase in a, phi, and phi^2
Xmat <- cbind(TT$a[usedObs], TT$phi[usedObs], TT$phi[usedObs]^2)
OneStdIncrease_x <- apply(Xmat, 2, sd, na.rm = TRUE) * coef(reg)[2:4]
OneStd_y <- sd(RPemp[usedObs], na.rm = TRUE)
OneStdIncrease_x / OneStd_y

# Table 6, Panel A, col 7,8
tab_model
tab_data


## -----------------------------------------------------------------------------
## Panel B: Volatility
## -----------------------------------------------------------------------------

Volmodel <- TTmodel$vol
Volemp   <- TT$VIX

## Panel B, col 1,2: Model regression
### ------------------------------------
reg <- lm(Volmodel ~ a, data = TT)
tab_model <- hac_robust(reg)
# Extract the logical index of used observations
usedObs <- complete.cases(TT$a, Volmodel)
# Extract the corresponding dates from TT
datesUsed <- TT$date[usedObs]
# Get first and last date used in the regression
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

### Panel B, col 1,2: Data regression
### ------------------------------------
reg <- lm(Volemp ~ a, data = TT)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$a, Volemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel B, col 1,2
tab_model
tab_data


### Panel B, col 3,4: Model regression
### ------------------------------------
reg <- lm(Volmodel ~ I(phi^2), data = TT)
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$phi, Volmodel)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

### Panel B, col 3,4: Data regression
### ------------------------------------
reg <- lm(Volemp ~ I(phi^2), data = TT)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$phi, Volemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel B, col 3,4
tab_model
tab_data


### Panel B, col 5,6: Model regression
### ------------------------------------
reg <- lm(Volmodel ~ a + I(phi^2), data = TT)
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, Volmodel)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

### Panel B, col 5,6: Data regression
### ------------------------------------
reg <- lm(Volemp ~ a + I(phi^2), data = TT)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, Volemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel B, col 5,6
tab_model
tab_data


### Panel B, col 7,8: Model regression
### ------------------------------------
reg <- lm(Volmodel ~ a + phi + I(phi^2), data = TT)
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, Volmodel)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)
OneStdIncrease_x <- apply(cbind(TT$a[usedObs], TT$phi[usedObs], TT$phi[usedObs]^2), 2, sd, na.rm = TRUE) *
  coef(reg)[2:4]
OneStd_y <- sd(Volmodel[usedObs], na.rm = TRUE)
OneStdIncrease_x / OneStd_y

### Panel B, col 7,8: Data regression
### ------------------------------------
reg <- lm(Volemp ~ a + phi + I(phi^2), data = TT)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, Volemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)
OneStdIncrease_x <- apply(cbind(TT$a[usedObs], TT$phi[usedObs], TT$phi[usedObs]^2), 2, sd, na.rm = TRUE) *
  coef(reg)[2:4]
OneStd_y <- sd(Volemp[usedObs], na.rm = TRUE)
OneStdIncrease_x / OneStd_y

# Table 6, Panel B, col 7,8
tab_model
tab_data


## =============================================================================
## 14. Table 7: Market Risk Premium and Return Volatility (X) vs. Model-Free Proxies (Y)
## =============================================================================

## -----------------------------------------------------------------------------
## Panel A: Market risk premium vs. model-free proxies
## -----------------------------------------------------------------------------
# φ_proxy is inflation deviations from its mean; 'na.rm = TRUE' ignores NaNs when computing the mean.
phi_proxy <- TT$inflation - mean(TT$inflation, na.rm = TRUE)

# Pearson correlation between TT$phi and phi_proxy; use only rows with no NaNs in either series.
cor(TT$phi, phi_proxy, use = "complete.obs")

# compute long-term inflation expectation: mean across horizons 5 to 30
ExpectedInflationMat <- NULL
for (i in 5:30) { # i = [5, 6, 7, ..., 30]
  varname <- sprintf("%d year Expected Inflation", i)
  ExpectedInflationMat <- cbind(ExpectedInflationMat, TT[[varname]]) # Horizontally concatenate each series into the matrix (aligned by TT's time index).
}
# Compute the mean of the expected inflation matrix across the specified horizons
LTExpectedInflationMean <- rowMeans(ExpectedInflationMat, na.rm = TRUE)
# 5th percentile of inflation as a "low-inflation threshold".
LowQuantileInflation <- quantile(TT$inflation, 0.05, na.rm = TRUE)
# Build a_proxy. Denominator uses elementwise max(low-quantile threshold, actual inflation) to avoid explosive ratios at very low inflation.
a_proxy <- -LTExpectedInflationMean / pmax(LowQuantileInflation, TT$inflation, na.rm = TRUE)
# Correlate model-implied a (TT$a) with the proxy a_proxy using complete cases.
cor(TT$a, a_proxy, use = "complete.obs")


### Table 7, Panel A, col. 1
### ------------------------------------

# OLS regression: RPemp_t = β0 + β1 * a_proxy_t + ε_t
reg <- lm(RPemp ~ a_proxy)
tab_data <- hac_robust(reg)   # HAC-robust (Newey–West) SEs
tab_data

# Identify used observations (exclude NA)
usedObs <- complete.cases(RPemp, a_proxy)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

### Table 7, Panel A, col. 2
### ------------------------------------
reg <- lm(RPemp ~ I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data

### Table 7, Panel A, col. 3
### ------------------------------------
reg <- lm(RPemp ~ a_proxy + I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data

### Table 7, Panel A, col. 4
### ------------------------------------
reg <- lm(RPemp ~ a_proxy + phi_proxy + I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data


## -----------------------------------------------------------------------------
## Panel B: Market return volatility vs. model-free proxies
## -----------------------------------------------------------------------------

### Table 7, Panel B, col. 1
### ------------------------------------
reg <- lm(Volemp ~ a_proxy)
tab_data <- hac_robust(reg)
tab_data

usedObs <- complete.cases(Volemp, a_proxy)
usedObsBench <- usedObs  # store benchmark logical index
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

### Table 7, Panel B, col. 2
### ------------------------------------
reg <- lm(Volemp[usedObsBench] ~ I(phi_proxy[usedObsBench]^2))
tab_data <- hac_robust(reg)
tab_data

### Table 7, Panel B, col. 3
### ------------------------------------
reg <- lm(Volemp ~ a_proxy + I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data

### Table 7, Panel B, col. 4
### ------------------------------------
reg <- lm(Volemp ~ a_proxy + phi_proxy + I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data


## =============================================================================
## 15. Table 8: Impact of Monetary Uncertainty on Market Risk Premium and Volatility
## =============================================================================

# nua_proxy = monetary policy uncertainty (Bloom index)
nua_proxy <- TT$MPU1   # monetary policy uncertainty, Bloom
IndTime <- !is.na(nua_proxy)   # IndTime marks the non-missing (valid) observations.

# Define time periods with high monetary policy uncertainty.
Times_HighnuaTrue <- as.numeric(TT$nua > quantile(TT$nua[!is.na(TT$nua)], 0.5, na.rm = TRUE))
Times_Highnua <- as.numeric(nua_proxy > quantile(nua_proxy[!is.na(nua_proxy)], 0.5, na.rm = TRUE))
Times_Highnua[is.na(nua_proxy)] <- NA  # Set NA where nua_proxy is NA

## -----------------------------------------------------------------------------
## Panel A: State Variables
## -----------------------------------------------------------------------------

### Table 8, Panel A, col. 1
### ------------------------------------
reg <- lm(RPemp ~ TT$a + TT$phi + I(TT$phi^2) + Times_HighnuaTrue)
tab_data <- hac_robust(reg)
tab_data

### Table 8, Panel A, col. 2
### ------------------------------------
reg <- lm(RPemp ~ TT$a + TT$phi + I(TT$phi^2) + Times_HighnuaTrue + I((TT$phi^2) * Times_HighnuaTrue))
tab_data <- hac_robust(reg)
tab_data

### Table 8, Panel A, col. 3 
### ------------------------------------
reg <- lm(Volemp ~ TT$a + TT$phi + I(TT$phi^2) + Times_HighnuaTrue)
tab_data <- hac_robust(reg)
tab_data

### Table 8, Panel A, col. 4
### ------------------------------------
reg <- lm(Volemp ~ TT$a + TT$phi + I(TT$phi^2) + Times_HighnuaTrue + I((TT$phi^2) * Times_HighnuaTrue))
tab_data <- hac_robust(reg)
tab_data


## -----------------------------------------------------------------------------
## Panel B: Model-Free Proxies
## -----------------------------------------------------------------------------

### Table 8, Panel B, col. 1
### ------------------------------------
reg <- lm(RPemp ~ a_proxy + phi_proxy + I(phi_proxy^2) + Times_Highnua)
tab_data <- hac_robust(reg)
tab_data

### Table 8, Panel B, col. 2
### ------------------------------------
reg <- lm(RPemp ~ a_proxy + phi_proxy + I(phi_proxy^2) + Times_Highnua + I((phi_proxy^2) * Times_Highnua))
tab_data <- hac_robust(reg)
tab_data

### Table 8, Panel B, col. 3
### ------------------------------------
reg <- lm(Volemp ~ a_proxy + phi_proxy + I(phi_proxy^2) + Times_Highnua)
tab_data <- hac_robust(reg)
tab_data

### Table 8, Panel B, col. 4
### ------------------------------------
reg <- lm(Volemp ~ a_proxy + phi_proxy + I(phi_proxy^2) + Times_Highnua + I((phi_proxy^2) * Times_Highnua))
tab_data <- hac_robust(reg)
tab_data


# ==============================================================================
# 16. Predictive regressions for future realized excess returns
# ==============================================================================

# set n equal to 1 or 5 or 10
n <- 1  # horizon (in years) in predictive regressions

# Compute annualized n-year-forward moving realized excess return
TT$ForwardRealizedRP <- 12 * zoo::rollmean(TT$exret, k = n * 12, align = "left", fill = NA)  # equivalent to movmean([0 n*12-1])
ForwardRealizedRP <- c(TT$ForwardRealizedRP[-1], NA)  # shift forward by 1, last value set to NA

## Table 9, Panel A or B or C (depending on choice of n above), col. 1
## --------------------------------------------------------------------
reg <- lm(ForwardRealizedRP ~ a_proxy)
tab_data <- hac_robust(reg)
tab_data

# Extract the logical index of used observations
usedObs <- complete.cases(a_proxy, ForwardRealizedRP)
usedObsBench <- usedObs

# Extract the corresponding dates from TT
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)


## Table 9, Panel A or B or C (depending on choice of n above), col. 2
## --------------------------------------------------------------------
reg <- lm(ForwardRealizedRP[usedObsBench] ~ I(phi_proxy[usedObsBench]^2))
tab_data <- hac_robust(reg)
tab_data


## Table 9, Panel A or B or C (depending on choice of n above), col. 3
## ---------------------------------------------------------------------
reg <- lm(ForwardRealizedRP ~ a_proxy + I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data


## Table 9, Panel A or B or C (depending on choice of n above), col. 4
## --------------------------------------------------------------------
reg <- lm(ForwardRealizedRP ~ a_proxy + phi_proxy + I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data