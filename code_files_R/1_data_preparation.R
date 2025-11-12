# Load packages
library(tidyverse)   # Core data manipulation & visualization (ggplot2, dplyr, tidyr, etc.)
library(dplyr)       # Data manipulation: filter, select, mutate, summarise, join
library(tidyr)       # Data reshaping (pivot_longer, pivot_wider, etc.)
library(readxl)      # Read Excel (.xlsx) files
library(readr)       # Read delimited text files (CSV, TSV, etc.)
library(zoo)         # Time series operations and rolling/moving functions
library(tsibble)     # Modern time series and time-indexed tibbles
library(lubridate)   # Date and time handling
library(xts)         # Extended time series structure for merging & alignment

# Clear
rm(list = ls())  # Clear all variables
graphics.off()   # Close all figures
cat("\014")      # Clear console

# Global constant
Delta <- 1/12  # data frequency (monthly)

# Path
setwd("/Users/yarusu/Library/CloudStorage/OneDrive-UCSanDiego/DS/ML project/AndreiHaslerJFEReplication_V2/MatlabCodeAndData/data")

# ==============================================================================
# Extract data and create timetable
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Output Gap
# ------------------------------------------------------------------------------
DataGap <- read_excel("OutputGapData.xlsx")
colnames(DataGap) <- c("date", "OutputGapQSimple")

# Convert to Date and shift +2 months (to middle month of quarter)
DataGap$date <- ymd(DataGap$date) %m+% months(2)

# Convert to continuously compounded form
DataGap <- DataGap %>%
  mutate(OutputGapQ = log(1 + OutputGapQSimple / 100))

# Define time range
Beg <- min(DataGap$date)
End <- ymd("2023-12-01")

# Create full monthly sequence
DataGapMonthly <- DataGap %>%
  complete(date = seq.Date(Beg, End, by = "month")) %>%
  fill(OutputGapQ, .direction = "up") %>%
  select(date, OutputGapQ) %>%
  filter(date <= ymd("2023-12-01"))

colnames(DataGapMonthly) <- c("date", "OutputGap")

# ------------------------------------------------------------------------------
# 2. Federal Funds Rate
# ------------------------------------------------------------------------------
DataFedFunds <- read_excel("FEDFUNDS.xlsx", skip = 7)
colnames(DataFedFunds) <- c("date", "fedfundsimple") # Rename variable
DataFedFunds$date <- as.Date(DataFedFunds$date) # Convert date column to Date format
DataFedFunds$fedfund <- log(1 + DataFedFunds$fedfundsimple / 100) # Compute continuously compounded (log) annualized rate

# ------------------------------------------------------------------------------
# 3. Goyal & Welch data
# ------------------------------------------------------------------------------
DataGoyal1 <- read_excel("GoyalWelch.xlsx")

# Extract year and month from the first column (YYYYMM format)
dateStr <- sprintf("%06d", DataGoyal1$yyyymm) # Convert numeric date to string
years <- as.numeric(substr(dateStr, 1, 4))
months <- as.numeric(substr(dateStr, 5, 6))
datetimevector <- as.Date(paste(years, months, "01", sep = "-")) # Create datetime vector

# Combine datetime with the remaining columns
DataGoyal <- data.frame(
  date = datetimevector,
  DataGoyal1[, -1]
)

# Assign variable names
colnames(DataGoyal) <- c(
  "date","index","d12","e12","bm","tbl","aaa","baa",
  "lty","ntis","rfree","infl","ltr","corpr","svar",
  "csp","sp_vw","sp_vwx"
)

# Compute 12-month simple moving average of log(1+infl) (annualized continuously compounded inflation)
DataGoyal$infl <- as.numeric(DataGoyal$infl)
DataGoyal$sp_vw <- as.numeric(DataGoyal$sp_vw)
DataGoyal <- DataGoyal %>%
  mutate(
    inflation = 12 * rollapply(log(1 + infl), 12, mean, fill = NA, align = "right"),
    retnom = log(1 + sp_vw)  # continuously compounded S&P500 return
  )

# Keep selected variables
DataGoyal <- DataGoyal %>%
  select(date, retnom, d12, index, inflation)

# ------------------------------------------------------------------------------
# 4. Real GDP data
# ------------------------------------------------------------------------------
DataGDP <- read_excel("RealGDP.xlsx")

# Rename columns
DataGDP <- DataGDP %>%
  rename(
    date = 1,   # first column is date/year
    gdp  = 2    # second column is GDP
  )

# Compute continuously compounded annual growth rate
DataGDP <- DataGDP %>%
  arrange(date) %>%  # ensure data is sorted by date
  mutate(
    GDPGrowthY = c(NA, log(gdp[-1] / gdp[-n()])),  # [NaN; log(gdp(2:end)/gdp(1:end-1))]
    date = as.Date(paste0(year(date), "-12-01"))   # assign date to December
  )

# Convert annual GDP growth to monthly
Beg <- min(DataGDP$date, na.rm = TRUE)
End <- max(DataGDP$date, na.rm = TRUE)
all_months <- seq(Beg, End, by = "month") # Generate monthly sequence

# Fill each month with the corresponding annual GDP growth
DataGDPMonthly <- data.frame(date = all_months) %>%
  left_join(DataGDP %>% select(date, GDPGrowthY), by = "date") %>%
  fill(GDPGrowthY, .direction = "up") %>%   
  rename(GDPGrowth = GDPGrowthY) 

# ==============================================================================
# Download additional variables
# ==============================================================================

# ------------------------------------------------------------------------------
# 5. Monthly Financial Risk and Equity Premium Data
# ------------------------------------------------------------------------------
DataMartinDaily <- read_csv("epbound_Martin.csv")  # daily
DataChabiYoDaily <- read_csv("epbound_ChabiYo.csv") # daily
DataVIXDaily <- read_csv("VIX.csv")                 # daily

# Ensure the date column is Date type
DataMartinDaily <- DataMartinDaily %>%
  mutate(date = mdy(date))  

DataChabiYoDaily <- DataChabiYoDaily %>%
  mutate(date = mdy(date))

DataVIXDaily <- read_csv("VIX.csv") %>%
  mutate(date = ymd(date))

# Convert to monthly (last value of the month)
# Martin ERP
DataRPMartinMonthly <- DataMartinDaily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(across(-date, ~ last(.x[!is.na(.x)]), .names = "{.col}"), .groups = "drop") %>%
  rename(date = month)

# Chabi-Yo ERP
DataRPChabiYoMonthly <- DataChabiYoDaily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(across(-date, ~ last(.x[!is.na(.x)]), .names = "{.col}"), .groups = "drop") %>%
  rename(date = month)

# VIX
DataVIXMonthly <- DataVIXDaily %>%
  mutate(month = floor_date(date, "month")) %>%     
  group_by(month) %>%
  summarise(VIX = last(VIX[!is.na(VIX)]), .groups = "drop") %>%  
  rename(date = month) %>%
  mutate(VIX = VIX / 100)

# ------------------------------------------------------------------------------
# 6. GDP growth forecast (median)
# ------------------------------------------------------------------------------
GDPgrowthForecastData <- read_csv("Median_RGDP_Growth.csv")

# Convert YEAR and QUARTER to a date (assuming quarters end in Mar, Jun, Sep, Dec)
GDPgrowthForecastData <- GDPgrowthForecastData %>%
  mutate(date = make_date(YEAR, QUARTER * 3, 1))

# Convert to a tibble with date as row identifier (similar to timetable)
DataGDPgrowthForecast <- GDPgrowthForecastData %>%
  select(-YEAR, -QUARTER) %>%  # Drop the YEAR and QUARTER columns (we already have 'date')
  column_to_rownames("date") %>%  # Set 'date' as the row names (like using time as index)
  as_tibble(rownames = "date") %>%  # Convert back to a tibble and keep 'date' as a column
  mutate(date = as.Date(date))

# Create monthly data: take the next available value for each month
# First, generate a full sequence of monthly dates
beg <- min(DataGDPgrowthForecast$date)
end <- max(DataGDPgrowthForecast$date)
monthly_dates <- seq.Date(beg, end, by = "month")

# Merge monthly dates with quarterly data, carry last observation forward
DataGDPgrowthForecastMonthly <- tibble(date = monthly_dates) %>%
  left_join(DataGDPgrowthForecast, by = "date") %>%
  arrange(date) %>%
  fill(DRGDP2, .direction = "down")  # last value carried forward

# Convert to continuously compounded growth
DataGDPgrowthForecastMonthly <- DataGDPgrowthForecastMonthly %>%
  mutate(mudelta = log(1 + DRGDP2 / 100)) %>%
  select(date, mudelta)

# ------------------------------------------------------------------------------
# 7. Inflation Expectations Cleveland Fed
# ------------------------------------------------------------------------------

# Dates at the end have a different format for whatever reason, hence the code below
# Detect import options and read CSV with first column as string
DataInflationExpectationMonthly <- read_csv(
  "InflationExpectationsClevelandFed.csv",
  col_types = cols(
    date = col_character()
  )
)

# First attempt: ISO format (e.g., 2022-09-01)
dates1 <- ymd(DataInflationExpectationMonthly$date, quiet = TRUE)

# Second attempt: short format (e.g., 1-1-2023) for failed cases
bad_idx <- is.na(dates1)
dates2 <- mdy(DataInflationExpectationMonthly$date[bad_idx], quiet = TRUE)

# Fill missing values
dates1[bad_idx] <- dates2

# Replace date column with parsed dates
DataInflationExpectationMonthly <- DataInflationExpectationMonthly %>%
  mutate(date = dates1) %>%
  arrange(date)  # ensure chronological order

# Convert the data frame to an xts time-series object
DataInflationExpectationMonthly$date <- as.Date(DataInflationExpectationMonthly$date)
DataInflationExpectationMonthly_xts <- xts(
  DataInflationExpectationMonthly %>% select(-date),
  order.by = DataInflationExpectationMonthly$date
)

# ------------------------------------------------------------------------------
# 8. Monetary Policy Uncertainty, MPU
# ------------------------------------------------------------------------------
MPUData <- read.csv("MPU.csv")

# Create date column (first day of each month)
MPUData$date <- as.Date(
  paste(MPUData$Year, MPUData$Month, "01", sep = "-"),
  format = "%Y-%m-%d"
)

# Convert to tibble and clean up
DataMPUMonthly <- MPUData %>%
  select(-Year, -Month) %>%   # remove Year and Month columns
  mutate(across(-date, ~ .x / 100)) %>%  # divide all numeric variables by 100
  drop_na() 

# ------------------------------------------------------------------------------
# 9. 5-Year CPI Forecast Data
# ------------------------------------------------------------------------------
CPI5YRForecastData <- read.csv("Median_CPI5YR.csv")

# Create date column assuming quarters end in Mar, Jun, Sep, Dec
CPI5YRForecastData <- CPI5YRForecastData %>%
  mutate(date = make_date(YEAR, QUARTER * 3, 1)) %>%
  select(-YEAR, -QUARTER) # Convert to tibble and clean up

# Define monthly sequence from first to last date
Beg <- min(CPI5YRForecastData$date, na.rm = TRUE)
End <- max(CPI5YRForecastData$date, na.rm = TRUE)

# Convert quarterly to monthly data — take next available value (same logic as MATLAB 'retime(..., "next")')
DataCPI5YRForecastMonthly <- CPI5YRForecastData %>%
  complete(date = seq(Beg, End, by = "month")) %>% # fill monthly sequence
  arrange(date) %>%
  fill(CPI5YR, .direction = "up") %>%   # propagate next available value
  mutate(CPI5YR = CPI5YR / 100)   # divide by 100 (continuously compounded style)

# ==============================================================================
# Merge All Monthly Datasets into One Unified DataFrame
# ==============================================================================

# List of all monthly datasets
TT_list <- list(
  DataGDPMonthly,
  DataFedFunds,
  DataGoyal,
  DataGapMonthly,
  DataVIXMonthly,
  DataRPChabiYoMonthly,
  DataRPMartinMonthly,
  DataGDPgrowthForecastMonthly,
  DataInflationExpectationMonthly,
  DataMPUMonthly,
  DataCPI5YRForecastMonthly
)


# Initialize with the first dataset
TTFull <- TT_list[[1]]

# Merge all datasets by date (full outer join)
for (i in 2:length(TT_list)) {
  TTFull <- full_join(TTFull, TT_list[[i]], by = "date")
}

# Restrict to the sample period (1954-07-01 to 2023-12-01)
TT <- TTFull %>%
  filter(date >= as.Date("1954-07-01") & date <= as.Date("2023-12-01")) %>%
  arrange(date)

# Compute other useful time series
TT <- TT %>%
  mutate(
    rfnom = fedfund,                        # nominal risk-free rate
    exret = retnom - rfnom / 12,            # excess stock return (S&P500)
    rfreal = rfnom - inflation,             # real risk-free rate
    dp = log(d12 / index),                  # log dividend-price ratio
    pd = -dp                                # log price-dividend ratio
  )

# ==============================================================================
# Descriptive statistics
# ==============================================================================

# ------------------------------------------------------------------------------
# 1.	Dataset Overview
# ------------------------------------------------------------------------------

# Sample Size
sample_info <- data.frame(
  Metric = c("Start Date", "End Date", "Number of Observations", "Number of Variables"),
  Value = c(
    as.character(min(TT$date, na.rm = TRUE)),
    as.character(max(TT$date, na.rm = TRUE)),
    nrow(TT),
    ncol(TT)
  )
)
sample_info

# Missingness Summary
library(knitr)

na_summary <- TT %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Missing_Count") %>%
  mutate(
    Missing_Percent = round(Missing_Count / nrow(TT) * 100, 2)
  ) %>%
  arrange(desc(Missing_Percent))

kable(na_summary, caption = "Missing Values for All 64 Variables in TT Dataset")

# ------------------------------------------------------------------------------
# 2.	Variable-Level Statistics
# ------------------------------------------------------------------------------
library(moments)

TT_numeric <- TT %>% select(where(is.numeric))

desc_stats <- lapply(TT_numeric, function(x) {
  data.frame(
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    skewness = ifelse(sum(!is.na(x)) > 2, skewness(x, na.rm = TRUE), NA),
    kurtosis = ifelse(sum(!is.na(x)) > 3, kurtosis(x, na.rm = TRUE), NA)
  )
}) %>%
  bind_rows(.id = "Variable")

desc_stats <- desc_stats %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

desc_stats
