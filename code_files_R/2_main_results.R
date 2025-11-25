library(tidyverse)
library(lubridate)
library(patchwork)
library(numDeriv)
library(stats)      
library(sandwich)  
library(lmtest)     
library(scales)
library(ggplot2)

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

# Fix Panel Order: GDP -> Output Gap -> Inflation -> Fed Funds
TT_long$variable <- factor(TT_long$variable, 
                           levels = c("GDPGrowth", "OutputGap", "inflation", "rfnom"))

facet_labels <- c(
  GDPGrowth = "GDP growth",
  OutputGap = "Output gap",
  inflation = "Inflation",
  rfnom = "Fed funds rate"
)

# Annotations
fed_chairs <- data.frame(
  Chair = c("Martin", "Burns", "Miller", "Volcker", "Greenspan", "Bernanke", "Yellen", "Powell"),
  Start = as.Date(c("1951-04-02", "1970-02-01", "1978-03-08", "1979-08-06", "1987-08-11", "2006-02-01", "2014-02-03", "2018-02-05")),
  End   = as.Date(c("1970-01-31", "1978-03-07", "1979-08-05", "1987-08-11", "2006-01-31", "2014-01-31", "2018-02-03", "2023-12-31"))
)

fed_chairs$Start <- pmax(fed_chairs$Start, min(TT$date))
fed_chairs$End   <- pmin(fed_chairs$End, max(TT$date))

shading_data <- fed_chairs[c(2, 4, 6, 8), ]

fed_chairs$Midpoint <- fed_chairs$Start + (fed_chairs$End - fed_chairs$Start)/2

chair_labels <- data.frame(
  date = fed_chairs$Midpoint,
  value = 0.08, 
  variable = factor("GDPGrowth", levels = levels(TT_long$variable)), 
  label = fed_chairs$Chair
)

# Plotting
p <- ggplot() +
  geom_rect(data = shading_data, 
            aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4) + 
  geom_line(data = TT_long, aes(x = date, y = value), 
            color = "blue", size = 0.6) +
  geom_text(data = chair_labels, aes(x = date, y = value, label = label),
            size = 2.5, color = "black") +
  facet_grid(variable ~ ., scales = "free_y", switch = "y",
             labeller = labeller(variable = facet_labels)) +
  scale_x_date(expand = c(0, 0), date_breaks = "10 years", date_labels = "%Y") +
  scale_y_continuous(n.breaks = 4) + 
  labs(x = NULL, y = NULL, 
       title = "Figure 1: GDP Growth, Output Gap, Inflation, and Federal Funds Rate") +
  theme_bw() +
  theme(
    strip.background = element_blank(), 
    strip.placement = "outside",        
    strip.text.y.left = element_text(size = 10, color = "black", angle = 90), 

    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 10)),

    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey"),

    panel.border = element_rect(color = "black", fill = NA)
  )

print(p)
ggsave("Figure1_Replication.png", p, width = 10, height = 8, dpi = 300)

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
# Requires: LikelihoodFunc_Rate, GetVarMatrixParam_Rate, rTaylor

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
# 5. Inflation and Transmission Mechanism: ML estimation of inflation parameters (sigmapi, pibar, lambdapi, sigmaa, lambdaa)
# ==============================================================================
# Requires: LikelihoodFunc, GetVarMatrixParam, TT (timetable), Delta, rfTaylor, rNbar_estimated

# ------------------------------------------------------------------------------
# A. Setup Initial Parameters and Bounds
# ------------------------------------------------------------------------------
param0 <- c(0.0131, 0.0338, 0.2436, 0.7434, 0.9365)
lower_bounds <- c(0.001, -0.05, 0.01, 0.01, 0.01)
upper_bounds <- c(0.100,  0.20, 2.00, 5.00, 5.00)

# ------------------------------------------------------------------------------
# B. Estimate Parameters via Maximum Likelihood
# ------------------------------------------------------------------------------
# Define wrapper function for optim (returns only the scalar objective value)
obj_func <- function(x) {
  LikelihoodFunc(TT$inflation, TT$fedfund, x)$sumloglik
}

cat("Starting Maximum Likelihood Estimation using L-BFGS-B...\n")

# Run optimization
result <- optim(
  par = param0,
  fn = obj_func,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(factr = 1e7, maxit = 2000) # High precision, allow more iterations
)

param <- result$par
fval <- result$value

# ------------------------------------------------------------------------------
# C. Post-Estimation Computations
# ------------------------------------------------------------------------------
# Recompute likelihood and filtered states using the final estimated parameters
out <- LikelihoodFunc(TT$inflation, TT$fedfund, param)
sumloglik <- out$sumloglik
logLik <- out$logLik
a <- out$a                # Filtered perceived effectiveness (a_t)
nua <- out$nua            # Filtered uncertainty (nu_t)
LTMpi <- out$LTMpi        # Long-term mean inflation
stdeps <- out$stdeps      # Standardized residuals

# Compute variance-covariance matrix of parameters (Hessian based)
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
cat("ML ESTIMATION RESULTS: INFLATION & LEARNING BLOCK\n")
cat("Parameter Vector: [sigmapi, pibar, lambdapi, sigmaa, lambdaa]\n")
cat("Estimated Parameters:\n")
print(param)
cat("\nStandard Errors:\n")
print(param_std)
cat("\nT-Statistics:\n")
print(param_tstats)
cat("\nP-Values:\n")
print(param_pvalues)

# ------------------------------------------------------------------------------
# D. Update Data and State Variables
# ------------------------------------------------------------------------------
# Insert ML-Kalman estimated state variables into the timetable
TT$a     <- a            # Perceived transmission effectiveness
TT$nua   <- nua          # Transmission uncertainty
TT$LTMpi <- LTMpi        # Perceived long-term inflation drift

# Estimated parameters from previous steps (Table 3 & Section 3.5)
beta_pi_est <- 0.39908112  # Taylor Rule response to inflation
beta_y_est  <- 0.11391485  # Taylor Rule response to output gap
pibar_est   <- 0.03381797  # Long-run mean inflation

# Compute phi_t = beta_pi * (pi_t - pi_bar) + beta_y * y_t
TT$phi <- beta_pi_est * (TT$inflation - pibar_est) + beta_y_est * TT$OutputGap

# ------------------------------------------------------------------------------
# E. Empirical Analysis of Estimated Variables
# ------------------------------------------------------------------------------
# (a) Relationship between changes in a_t and inflation
change_a <- diff(TT$a)
change_pi <- diff(TT$inflation)
dummy <- as.numeric(TT$inflation > mean(TT$inflation)) # High inflation regime dummy

# Correlation in high inflation periods (Expected: Negative)
corr_high <- cor(change_a[dummy[-1] == 1], change_pi[dummy[-1] == 1])
# Correlation in low inflation periods (Expected: Positive)
corr_low <- cor(change_a[dummy[-1] == 0], change_pi[dummy[-1] == 0])

# (b) Moments of phi (monetary policy deviation)
phi_lead <- TT$phi[-1]
phi_lag <- TT$phi[-length(TT$phi)]
reg_autocorr_phi <- lm(phi_lead ~ phi_lag)
autocorrphi <- coef(reg_autocorr_phi)[2]
Volphi <- sd(TT$phi, na.rm = TRUE)

# (c) Moments of long-term mean of inflation (LTMpi)
LTMpi_lead <- LTMpi[-1]
LTMpi_lag <- LTMpi[-length(LTMpi)]
reg_LTMpi <- lm(LTMpi_lead ~ LTMpi_lag)
autocorLTMpi <- coef(reg_LTMpi)[2]
volLTMpi <- sd(LTMpi, na.rm = TRUE)

# (d) Correlation with 5-year inflation expectations
corr_CPI5YR <- cor(TT$CPI5YR, TT$LTMpi, use = "complete.obs")
reg_CPI5YR <- lm(CPI5YR ~ LTMpi, data = TT)

# HAC robust regression (Newey-West)
library(sandwich)
library(lmtest)
hac_results <- coeftest(reg_CPI5YR, vcov = NeweyWest(reg_CPI5YR, prewhite = FALSE))

# Display analysis results
cat("DIAGNOSTIC STATISTICS\n")
cat("Correlation (Delta_a, Delta_pi | High Inflation):", corr_high, "\n")
cat("Correlation (Delta_a, Delta_pi | Low Inflation): ", corr_low, "\n")
cat("\nPolicy Stance (phi):\n")
cat("  Autocorrelation:", autocorrphi, "\n")
cat("  Std. Deviation: ", Volphi, "\n")
cat("\nLong-Term Inflation Mean (LTMpi):\n")
cat("  Autocorrelation:", autocorLTMpi, "\n")
cat("  Std. Deviation: ", volLTMpi, "\n")
cat("\nValidation against Survey Data:\n")
cat("  Correlation(LTMpi, CPI5YR):", corr_CPI5YR, "\n")
cat("\nHAC-robust regression results (CPI5YR ~ LTMpi):\n")
print(hac_results)


# ==============================================================================
# 6. Figure 2: Investors' Estimate of the Transmission Coefficient
# ==============================================================================

# Define standard parameters, Order: sigmapi, pibar, lambdapi, sigmaa, lambdaa
my_params <- c(0.01290653, 0.03381797, 0.23475856, 0.71574256, 1.08948551)

# Rerun the Kalman Filter to generate the corrected a_t
out_paper <- LikelihoodFunc(TT$inflation, TT$fedfund, my_params)

# Update state variables in TT
TT$a <- out_paper$a
TT$LTMpi <- out_paper$LTMpifinal

# Prepare Fed Chair tenure data (for background shading)
# Define start and end dates for each Chair
fed_chairs <- data.frame(
  Chair = c("Martin", "Burns", "Miller", "Volcker", "Greenspan", "Bernanke", "Yellen", "Powell"),
  Start = as.Date(c("1951-04-02", "1970-02-01", "1978-03-08", "1979-08-06", "1987-08-11", "2006-02-01", "2014-02-03", "2018-02-05")),
  End   = as.Date(c("1970-01-31", "1978-03-07", "1979-08-05", "1987-08-11", "2006-01-31", "2014-01-31", "2018-02-03", "2023-12-31")) # End of sample
)

# Keep only Chairs within the sample period (post-1954)
fed_chairs$Start <- pmax(fed_chairs$Start, min(TT$date))
fed_chairs$End   <- pmin(fed_chairs$End, max(TT$date))

# Select every other Chair (2, 4, 6, 8) for alternating shading (Burns, Volcker, Bernanke, Powell)
# Or Martin, Miller, Greenspan, Yellen (light colors)
# In the paper's figure: Burns, Volcker, Bernanke, Powell have dark backgrounds
shading_data <- fed_chairs[c(2, 4, 6, 8), ]

# Set dual-axis scaling coefficient
coeff <- 6 

# Plotting
p <- ggplot() +
  # Plot background shading (Fed Chairs)
  geom_rect(data = shading_data, 
            aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.5) +
  
  # Plot left axis data: Transmission coefficient a_t (solid line, blue)
  geom_line(data = TT, aes(x = date, y = a), 
            color = "blue", size = 0.8) +
  
  # Plot right axis data: Inflation (dashed line, red)
  # Note: Multiply y by coeff to scale up
  geom_line(data = TT, aes(x = date, y = inflation * coeff), 
            color = "red", linetype = "dashed", size = 0.8) +
  
  # Set axes
  scale_y_continuous(
    # Left axis settings
    name = expression(paste("Transmission coefficient ", hat(a)[t])),
    limits = c(-1.1, 0.7), # Force limits to match the paper
    breaks = seq(-1.0, 0.5, by = 0.25),
    
    # Right axis settings (sec_axis)
    sec.axis = sec_axis(~ . / coeff, # Divide data back to display correct labels
                        name = expression(paste("Inflation ", pi[t])),
                        breaks = seq(-0.05, 0.15, by = 0.05))
  ) +
  
  # Set X axis
  scale_x_date(date_breaks = "10 years", date_labels = "%Y", 
               limits = c(as.Date("1954-01-01"), as.Date("2024-01-01")),
               expand = c(0,0)) +
  
  # Theme and aesthetics
  theme_bw() + 
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.y.left = element_text(color = "blue", size = 12, margin = margin(r = 10)),
    axis.text.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "red", size = 12, angle = 90, margin = margin(l = 10)),
    axis.text.y.right = element_text(color = "red"),
    axis.title.x = element_blank() # Remove Date label
  ) +
  
  # Add Chair names (Optional, manually adjust y position)
  annotate("text", x = as.Date("1960-01-01"), y = 0.65, label = "Martin", size = 3, color="grey40") +
  annotate("text", x = as.Date("1974-01-01"), y = 0.65, label = "Burns", size = 3, color="grey40") +
  annotate("text", x = as.Date("1978-09-01"), y = 0.65, label = "Miller", size = 3, color="grey40", angle=90) +
  annotate("text", x = as.Date("1983-01-01"), y = 0.65, label = "Volcker", size = 3, color="grey40") +
  annotate("text", x = as.Date("1996-01-01"), y = 0.65, label = "Greenspan", size = 3, color="grey40") +
  annotate("text", x = as.Date("2010-01-01"), y = 0.65, label = "Bernanke", size = 3, color="grey40") +
  annotate("text", x = as.Date("2016-01-01"), y = 0.65, label = "Yellen", size = 3, color="grey40") +
  annotate("text", x = as.Date("2021-01-01"), y = 0.65, label = "Powell", size = 3, color="grey40")

print(p)
ggsave("Figure2_Replication.png", p, width = 10, height = 5, dpi = 300)


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

# ------------------------------------------------------------------------------
# A. Data Summary Statistics
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# B. Model-Implied Time Series (from Mathematica output)
# ------------------------------------------------------------------------------
# Load model-implied data (from Mathematica export)
DataModel <- read.csv("data/ModelImpliedTimeSeriesNov2024.csv")

# Align model-implied data with empirical sample dates
DataModel$date <- TT$date   # ensure same timeline as empirical TT

# Convert to tibble for easy manipulation
library(tibbletime)
TTmodel <- as_tibble(DataModel)

# Compute model-implied moments (annualized, consistent with part (a))
meanrfreal_model <- mean(TTmodel$rfreal, na.rm = TRUE)
meanrfnom_model  <- mean(TTmodel$rfnom,  na.rm = TRUE)
meanSP_model     <- mean(TTmodel$RP,     na.rm = TRUE)
volSP_model      <- mean(TTmodel$vol,    na.rm = TRUE)

# Combine Data vs Model results into a summary table
Table4 <- data.frame(
  Moment = c("Real interest rate",
             "Nominal interest rate",
             "Market risk premium",
             "Market return volatility"),
  Data = c(meanrfreal, meanrfnom, meanSP, volSP),
  Model = c(meanrfreal_model, meanrfnom_model, meanSP_model, volSP_model)
)

# Display final table (rounded)
cat("\nTable 4: Asset-pricing Moments\n")
Table4[, 2:3] <- round(Table4[, 2:3], 4)
print(Table4)


# ==============================================================================
# 12. Table 5: Expected Output Growth, Real Interest Rate, and Log Price-Dividend Ratio (X) vs. Inflation and the Output Gap (Y).
# ==============================================================================

# ------------------------------------------------------------------------------
# (a) Col 1,2: Expected Output Growth
# ------------------------------------------------------------------------------
# Model regression
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

# Data regression
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

# ------------------------------------------------------------------------------
# (b) Col 3,4: Real Risk-Free Rate
# ------------------------------------------------------------------------------
# Model regression
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

# Data regression
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

# ------------------------------------------------------------------------------
# (c) Col 5,6: Log PD Ratio
# ------------------------------------------------------------------------------
# Model regression
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

# Data regression
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

# ------------------------------------------------------------------------------
# Panel A
# ------------------------------------------------------------------------------

# ------------------------------------
# Panel A, col 1,2: Model regression
# ------------------------------------
# corr between Martin(model) and ChabiYo RPs(data)
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

# ------------------------------------
# Panel A, col 1,2: Data regression
# ------------------------------------
reg <- lm(RPemp ~ TT$a)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$a, RPemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel A, col 1,2
tab_model
tab_data

# ------------------------------------
# Panel A, col 3,4: Model regression
# ------------------------------------
reg <- lm(RPmodel ~ I(TT$phi^2))
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$phi, RPmodel)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# ------------------------------------
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

# ------------------------------------
# Panel A, col 5,6: Model regression
# ------------------------------------
reg <- lm(RPmodel ~ TT$a + I(TT$phi^2))
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, RPmodel)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# ------------------------------------
# Panel A, col 5,6: Data regression
# ------------------------------------
reg <- lm(RPemp ~ TT$a + I(TT$phi^2))
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, RPemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel A, col 5,6
tab_model
tab_data

# ------------------------------------
# Panel A, col 7,8: Model regression
# ------------------------------------
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

# ------------------------------------
# Panel A, col 7,8: Data regression
# ------------------------------------
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

# ------------------------------------------------------------------------------
# Panel B: Volatility
# ------------------------------------------------------------------------------
Volmodel <- TTmodel$vol
Volemp   <- TT$VIX

# ------------------------------------
# Panel B, col 1,2: Model regression
# ------------------------------------
reg <- lm(Volmodel ~ a, data = TT)
tab_model <- hac_robust(reg)

# Extract the logical index of used observations
usedObs <- complete.cases(TT$a, Volmodel)

# Extract the corresponding dates from TT
datesUsed <- TT$date[usedObs]

# Get first and last date used in the regression
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# ------------------------------------
# Panel B, col 1,2: Data regression
# ------------------------------------
reg <- lm(Volemp ~ a, data = TT)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$a, Volemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel B, col 1,2
tab_model
tab_data

# ------------------------------------
# Panel B, col 3,4: Model regression
# ------------------------------------
reg <- lm(Volmodel ~ I(phi^2), data = TT)
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$phi, Volmodel)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# ------------------------------------
# Panel B, col 3,4: Data regression
# ------------------------------------
reg <- lm(Volemp ~ I(phi^2), data = TT)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$phi, Volemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel B, col 3,4
tab_model
tab_data

# ------------------------------------
# Panel B, col 5,6: Model regression
# ------------------------------------
reg <- lm(Volmodel ~ a + I(phi^2), data = TT)
tab_model <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, Volmodel)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# ------------------------------------
# Panel B, col 5,6: Data regression
# ------------------------------------
reg <- lm(Volemp ~ a + I(phi^2), data = TT)
tab_data <- hac_robust(reg)
usedObs <- complete.cases(TT$a, TT$phi, Volemp)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# Table 6, Panel B, col 5,6
tab_model
tab_data

# ------------------------------------
# Panel B, col 7,8: Model regression
# ------------------------------------
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

# ------------------------------------
# Panel B, col 7,8: Data regression
# ------------------------------------
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


# ==============================================================================
# 14. Table 7: Market Risk Premium and Return Volatility (X) vs. Model-Free Proxies (Y)
# ==============================================================================

# ------------------------------------------------------------------------------
# Panel A: Market risk premium vs. model-free proxies
# ------------------------------------------------------------------------------
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

# ------------------------------------
# Table 7, Panel A, col. 1
# ------------------------------------

# OLS regression: RPemp_t = β0 + β1 * a_proxy_t + ε_t
reg <- lm(RPemp ~ a_proxy)
tab_data <- hac_robust(reg)   # HAC-robust (Newey–West) SEs
tab_data

# Identify used observations (exclude NA)
usedObs <- complete.cases(RPemp, a_proxy)
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# ------------------------------------
# Table 7, Panel A, col. 2
# ------------------------------------
reg <- lm(RPemp ~ I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 7, Panel A, col. 3
# ------------------------------------
reg <- lm(RPemp ~ a_proxy + I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 7, Panel A, col. 4
# ------------------------------------
reg <- lm(RPemp ~ a_proxy + phi_proxy + I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------------------------------------------------
# Panel B: Market return volatility vs. model-free proxies
# ------------------------------------------------------------------------------

# ------------------------------------
# Table 7, Panel B, col. 1
# ------------------------------------
reg <- lm(Volemp ~ a_proxy)
tab_data <- hac_robust(reg)
tab_data

usedObs <- complete.cases(Volemp, a_proxy)
usedObsBench <- usedObs  # store benchmark logical index
datesUsed <- TT$date[usedObs]
firstDate <- min(datesUsed, na.rm = TRUE)
lastDate  <- max(datesUsed, na.rm = TRUE)

# ------------------------------------
# Table 7, Panel B, col. 2
# ------------------------------------
reg <- lm(Volemp[usedObsBench] ~ I(phi_proxy[usedObsBench]^2))
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 7, Panel B, col. 3
# ------------------------------------
reg <- lm(Volemp ~ a_proxy + I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 7, Panel B, col. 4
# ------------------------------------
reg <- lm(Volemp ~ a_proxy + phi_proxy + I(phi_proxy^2))
tab_data <- hac_robust(reg)
tab_data


# ==============================================================================
# 15. Table 8: Impact of Monetary Uncertainty on Market Risk Premium and Volatility
# ==============================================================================

# nua_proxy = monetary policy uncertainty (Bloom index)
nua_proxy <- TT$MPU1   # monetary policy uncertainty, Bloom
IndTime <- !is.na(nua_proxy)   # IndTime marks the non-missing (valid) observations.

# Define time periods with high monetary policy uncertainty.
Times_HighnuaTrue <- as.numeric(TT$nua > quantile(TT$nua[!is.na(TT$nua)], 0.5, na.rm = TRUE))
Times_Highnua <- as.numeric(nua_proxy > quantile(nua_proxy[!is.na(nua_proxy)], 0.5, na.rm = TRUE))
Times_Highnua[is.na(nua_proxy)] <- NA  # Set NA where nua_proxy is NA

# ------------------------------------------------------------------------------
# Panel A: State Variables
# ------------------------------------------------------------------------------

# ------------------------------------
# Table 8, Panel A, col. 1
# ------------------------------------
reg <- lm(RPemp ~ TT$a + TT$phi + I(TT$phi^2) + Times_HighnuaTrue)
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 8, Panel A, col. 2
# ------------------------------------
reg <- lm(RPemp ~ TT$a + TT$phi + I(TT$phi^2) + Times_HighnuaTrue + I((TT$phi^2) * Times_HighnuaTrue))
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 8, Panel A, col. 3 
# ------------------------------------
reg <- lm(Volemp ~ TT$a + TT$phi + I(TT$phi^2) + Times_HighnuaTrue)
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 8, Panel A, col. 4
# ------------------------------------
reg <- lm(Volemp ~ TT$a + TT$phi + I(TT$phi^2) + Times_HighnuaTrue + I((TT$phi^2) * Times_HighnuaTrue))
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------------------------------------------------
# Panel B: Model-Free Proxies
# ------------------------------------------------------------------------------

# ------------------------------------
# Table 8, Panel B, col. 1
# ------------------------------------
reg <- lm(RPemp ~ a_proxy + phi_proxy + I(phi_proxy^2) + Times_Highnua)
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 8, Panel B, col. 2
# ------------------------------------
reg <- lm(RPemp ~ a_proxy + phi_proxy + I(phi_proxy^2) + Times_Highnua + I((phi_proxy^2) * Times_Highnua))
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 8, Panel B, col. 3
# ------------------------------------
reg <- lm(Volemp ~ a_proxy + phi_proxy + I(phi_proxy^2) + Times_Highnua)
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 8, Panel B, col. 4
# ------------------------------------
reg <- lm(Volemp ~ a_proxy + phi_proxy + I(phi_proxy^2) + Times_Highnua + I((phi_proxy^2) * Times_Highnua))
tab_data <- hac_robust(reg)
tab_data


# ==============================================================================
# 16. Table 9: Predictive regressions for future realized excess returns
# ==============================================================================
# Set horizon (n = 1, 5, or 10)
n <- 10  

# Subset data to match sample period (Jan 1982 - Dec 2023)
idx_sample <- which(TT$date >= as.Date("1982-01-01") & TT$date <= as.Date("2023-12-31"))
TT_sub <- TT[idx_sample, ]

# Extract proxies for the specific sample
a_sub   <- a_proxy[idx_sample]
phi_sub <- phi_proxy[idx_sample]

# Compute annualized n-year-forward moving realized excess return
# 'partial = TRUE' keeps tail observations (simulates MATLAB movmean), ensuring N=503
window <- n * 12
raw_fwd <- zoo::rollapply(TT_sub$exret, width = window, FUN = mean, align = "left", partial = TRUE)

# Shift Y variable: predict return starting at t+1 using state at t
# Remove first element (t) and append NA at end
ForwardRealizedRP <- c(12 * raw_fwd[-1], NA)

# Create dataframe for regression (automatically handles alignment)
df <- data.frame(Y = ForwardRealizedRP, A = a_sub, Phi = phi_sub)
df <- na.omit(df)  # Remove the final NA row -> Obs = 503

# ------------------------------------
# Table 9, Panel A/B/C, col. 1
# ------------------------------------
reg <- lm(Y ~ A, data = df)
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 9, Panel A/B/C, col. 2
# ------------------------------------
reg <- lm(Y ~ I(Phi^2), data = df)
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 9, Panel A/B/C, col. 3
# ------------------------------------
reg <- lm(Y ~ A + I(Phi^2), data = df)
tab_data <- hac_robust(reg)
tab_data

# ------------------------------------
# Table 9, Panel A/B/C, col. 4
# ------------------------------------
reg <- lm(Y ~ A + Phi + I(Phi^2), data = df)
tab_data <- hac_robust(reg)
tab_data