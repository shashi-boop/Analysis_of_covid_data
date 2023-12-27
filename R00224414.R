
#SETTING THE WORKING DIRECTORY
setwd("C:\\Users\\shash\\OneDrive\\Desktop\\Time series and Factor Analysis")



#lOADING THE REQUIRED LIBRARY PACKAGES
library(tidyverse)
library(lubridate)
library(forecast)
library(ggplot2)
library(readr)
library(zoo)
#READING THE COVID DATA
data <- read.csv("C:\\Users\\shash\\OneDrive\\Desktop\\Time series and Factor Analysis\\COVID_DATASET.csv")
str(data)
# Filter for Chile and select relevant columns
chile_data <- data %>%
  filter(location == "Chile") %>%
  select(date, new_tests)

# Convert date column to a date object
chile_data$date <- as.Date(chile_data$date)

############################### DATA CLEANING AND MANIPULATION #################
# Interpolate missing values
str(chile_data)
chile_data <- chile_data %>%
  mutate(new_tests = approx(new_tests, n = n(), method = "linear", rule = 2)$y) %>%
  fill(new_tests)  
str(chile_data)
# Fill any remaining NA's at the beginning of the series
# Create weekly and monthly time series
chile_weekly <- chile_data %>%
  mutate(week = floor_date(date, "week")) %>%
  group_by(week) %>%
  summarize(weekly_tests = sum(new_tests, na.rm = TRUE))

chile_monthly <- chile_data %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarize(monthly_tests = sum(new_tests, na.rm = TRUE))

chile_daily <- chile_data %>%
  select(date, daily_tests = new_tests)
# Convert data frames to time series objects
chile_daily_ts <- ts(chile_daily$daily_tests, start = decimal_date(min(chile_daily$date)), frequency = 365.25)
chile_weekly_ts <- ts(chile_weekly$weekly_tests, start = decimal_date(min(chile_weekly$week)), frequency = 365.25 / 7)
chile_monthly_ts <- ts(chile_monthly$monthly_tests, start = decimal_date(min(chile_monthly$month)), frequency = 12)

# Plot daily, weekly, and monthly time series using autoplot
autoplot(chile_daily_ts) +
  ggtitle("Daily Tests Performed in Chile") +
  xlab("Date") + ylab("Number of Tests") +
  theme_minimal()+
  theme_bw()

autoplot(chile_weekly_ts) +
  ggtitle("Weekly Tests Performed in Chile") +
  xlab("Date") + ylab("Number of Tests") +
  theme_minimal()+
  theme_bw()

autoplot(chile_monthly_ts) +
  ggtitle("Monthly Tests Performed in Chile") +
  xlab("Date") + ylab("Number of Tests") +
  theme_minimal()

# Summary statistics
summary(chile_daily_ts)
summary(chile_weekly_ts)
################################# Preliminary descriptive analysis ############################
# Summary statistics
summary(chile_daily_ts)
summary(chile_weekly_ts)
summary(chile_monthly_ts)

autoplot(daily_ts) + ggtitle("Daily Tests Performed") + xlab("Time") + ylab("Number of Tests")
autoplot(weekly_ts) + ggtitle("Weekly Tests Performed") + xlab("Time") + ylab("Number of Tests")
autoplot(monthly_ts) + ggtitle("Monthly Tests Performed") + xlab("Time") + ylab("Number of Tests")
# Autocorrelation plots
acf(chile_daily_ts, main = "Daily Autocorrelation")
acf(chile_weekly_ts, main = "Weekly Autocorrelation")
acf(chile_monthly_ts, main = "Monthly Autocorrelation")
# Decompose daily, weekly, and monthly time series
chile_daily_decomp <- stl(chile_daily_ts, s.window = "periodic", robust = TRUE)
chile_weekly_decomp <- stl(chile_weekly_ts, s.window = "periodic", robust = TRUE)
chile_monthly_decomp <- stl(chile_monthly_ts, s.window = "periodic", robust = TRUE)

# Plot decompositions
autoplot(chile_daily_decomp) + ggtitle("Daily Tests Decomposition")
autoplot(chile_weekly_decomp) + ggtitle("Weekly Tests Decomposition")
autoplot(chile_monthly_decomp) + ggtitle("Monthly Tests Decomposition")

# Subset the time series by wave (replace the dates with the specific dates of each wave)
wave1 <- window(chile_daily_ts, start = decimal_date(as.Date("2020-03-15")), end = decimal_date(as.Date("2020-06-30")))
wave2 <- window(chile_daily_ts, start = decimal_date(as.Date("2020-07-01")), end = decimal_date(as.Date("2021-01-31")))
wave3 <- window(chile_daily_ts, start = decimal_date(as.Date("2021-02-01")), end = decimal_date(as.Date("2021-09-30")))
# Summary statistics
summary(wave1)
summary(wave2)
summary(wave3)

# Autocorrelation plots
acf(wave1, main = "Wave 1 Autocorrelation")
acf(wave2, main = "Wave 2 Autocorrelation")
acf(wave3, main = "Wave 3 Autocorrelation")

# Decompose time series for each wave
wave1_decomp <- stl(wave1, s.window = "periodic", robust = TRUE)
wave2_decomp <- stl(wave2, s.window = "periodic", robust = TRUE)
wave3_decomp <- stl(wave3, s.window = "periodic", robust = TRUE)

# Plot decompositions
autoplot(wave1_decomp) + ggtitle("Wave 1 Tests Decomposition")
autoplot(wave2_decomp) + ggtitle("Wave 2 Tests Decomposition")
autoplot(wave3_decomp) + ggtitle("Wave 3 Tests Decomposition")

################################## Time series Modelling #############################
library(forecast)

# Simple Exponential Smoothing (SES)
ses_model <- ses(chile_daily_ts)

# Holt's Linear Trend Method
holt_model <- holt(chile_daily_ts)

# TBATS model
tbats_model <- tbats(chile_daily_ts)

# Holt-Winters' Seasonal Method
hw_additive <- hw(chile_daily_ts, seasonal = "additive")
hw_multiplicative <- hw(chile_daily_ts, seasonal = "multiplicative")

# Compare the models using accuracy measures
accuracy(ses_model)
accuracy(holt_model)
accuracy(hw_additive)
accuracy(hw_multiplicative)
accuracy(tbats_model)

# Augmented Dickey-Fuller Test
library(tseries)
adf_test <- adf.test(chile_daily_ts)

# KPSS Test
kpss_test <- kpss.test(chile_daily_ts)

# Plot the time series
autoplot(chile_daily_ts)
# Differencing the time series
diff_ts <- diff(chile_daily_ts)

# Re-run the stationarity tests and plot the differenced series
adf_test_diff <- adf.test(diff_ts)
kpss_test_diff <- kpss.test(diff_ts)
autoplot(diff_ts)
# ACF and PACF plots
acf(diff_ts)
pacf(diff_ts)
# Fit the models
# Let auto.arima find the best model
model1 <- auto.arima(chile_daily_ts)

# Check the model summary
summary(model1)
model1 <- auto.arima(chile_daily_ts, order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = S))
model2 <- auto.arima(chile_daily_ts, order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = S))
model3 <- auto.arima(chle_daily_ts, order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = S))

# Residual analysis
autoplot(residuals(model1)) + ggtitle("Model 1 Residuals")
autoplot(residuals(model2)) + ggtitle("Model 2 Residuals")
autoplot(residuals(model3)) + ggtitle("Model 3 Residuals")


##################################### FORECAST ######################################
# Load forecast library
library(forecast)

# Fit Holt-Winters (Holt's) model with additive seasonality
hw_additive <- hw(chile_daily_ts, seasonal = "additive", h = 28)

# Fit Holt-Winters (Holt's) model with multiplicative seasonality
hw_multiplicative <- hw(chile_daily_ts, seasonal = "multiplicative", h = 28)

# Display the forecasts
plot(hw_additive, main = "Additive Holt-Winters Forecast")
plot(hw_multiplicative, main = "Multiplicative Holt-Winters Forecast")

# Fit auto.arima model
auto_model <- auto.arima(chile_daily_ts)

# Forecast 4 weeks (28 days) ahead
auto_forecast <- forecast(auto_model, h = 28)

# Display the forecast
plot(auto_forecast, main = "Auto ARIMA Forecast")
# Fit chosen ARIMA model
chosen_model <- Arima(chile_daily_ts, order = c(1, 1, 1))

# Forecast 4 weeks (28 days) ahead
chosen_forecast <- forecast(chosen_model, h = 28)

# Display the forecast
plot(chosen_forecast, main = "Chosen ARIMA Forecast")
# Remove last 28 days
chile_daily_ts_holdout <- window(chile_daily_ts, end = length(chile_daily_ts) - 28)

# Fit models with the holdout sample
hw_additive_holdout <- hw(chile_daily_ts_holdout, seasonal = "additive", h = 28)
hw_multiplicative_holdout <- hw(chile_daily_ts_holdout, seasonal = "multiplicative", h = 28)
auto_model_holdout <- auto.arima(chile_daily_ts_holdout)
chosen_model_holdout <- Arima(chile_daily_ts_holdout, order = c(1, 1, 1))

# Forecast 28 days ahead
hw_additive_forecast_holdout <- forecast(hw_additive_holdout, h = 28)
hw_multiplicative_forecast_holdout <- forecast(hw_multiplicative_holdout, h = 28)
auto_forecast_holdout <- forecast(auto_model_holdout, h = 28)
chosen_forecast_holdout <- forecast(chosen_model_holdout, h = 28)

# Compare forecasts with actual values
actual_values <- window(chile_daily_ts, start = length(chile_daily_ts) - 27)

# Calculate accuracy metrics
accuracy(hw_additive_forecast_holdout, actual_values)
accuracy(hw_multiplicative_forecast_holdout, actual_values)
accuracy(auto_forecast_holdout, actual_values)
accuracy(chosen_forecast_holdout, actual_values)
