
# Registration: 202454202008  -> last 3 digits = 008 -> studentID = 8
#                              -> last digit   = 8   -> Q2(8)


pkgs<-c('forecast','tseries','waveslim','rugarch','FinTS','lmtest')
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

install.packages(
  "D:\\M-1\\M-1-20260728T144156Z-1-001\\M-1\\Time series R code\\TS2025_0.1.0 2.tar.gz",
  repos=NULL,
  type="source"
)
library(TS2025)

ls("package:TS2025", all.names = TRUE)

sid<-54202008

# temp_ts function (example)
temp_ts=function(regID){
  library(forecast)
  # Set seed for reproducibility
  set.seed(sid)
  
  # Define the time series length: monthly data from Oct 2010 to Sep 2025
  start_year <- 2010
  start_month <- 10
  end_year <- 2025
  end_month <- 9
  
  n_years <- end_year - start_year + 1
  n_months <- n_years * 12 - (start_month - 1) - (12 - end_month)
  
  sim_data <- arima.sim(
    n = n_months,
    model = list(order = c(1,1,2), ar = 0.5, ma = c(0.3, -0.2))
  )
  
  round(ts(sim_data, start = c(start_year, start_month), frequency = 12))
}

#Helper fun:
to_ts <- function(x, start = NULL, frequency = 12) {
  if (is.ts(x)) return(x)
  if (is.data.frame(x)) {
    num <- x[sapply(x, is.numeric)]
    v <- num[[ncol(num)]]          # last numeric column = the series
  } else v <- as.numeric(x)
  if (is.null(start)) ts(v, frequency = frequency) else ts(v, start = start, frequency = frequency)
}




#question-1:

raw  <- temp_ts(sid)
str(raw); head(raw)
y    <- to_ts(raw)        # if the data has dates, add start = c(YYYY, M)
n    <- length(y); h <- 12
train <- window(y, end   = time(y)[n - h])
test  <- window(y, start = time(y)[n - h + 1])

# 1. Exploration
plot(y, main = "Monthly temperature", ylab = "Temp")
monthplot(y); seasonplot(y, year.labels = FALSE, main = "Seasonal plot")
plot(stl(train, s.window = "periodic"))
summary(train)

# 2. Stationarity
adf.test(train); kpss.test(train)
ndiffs(train); nsdiffs(train)
tsdisplay(train, main = "Training data")
tsdisplay(diff(train, 12))

# 3. Model fitting (seasonal ARIMA)
fit1 <- auto.arima(train, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
summary(fit1)
checkresiduals(fit1)                         # includes Ljung-Box
Box.test(residuals(fit1), lag = 24, type = "Ljung-Box",
         fitdf = length(coef(fit1)))
shapiro.test(residuals(fit1))

# 4. Forecast and accuracy
fc1 <- forecast(fit1, h = h)
plot(fc1); lines(test, col = "red")
acc1 <- accuracy(fc1, test)                  # ME, RMSE, MAE, MPE, MAPE, MASE
print(acc1)


# Check actual test values
print(test)

# Count zero values
sum(test == 0, na.rm = TRUE)

# Identify zero-value positions
which(test == 0)

# QUESTION 2 : Wavelet + Box-Jenkins (ARIMA) vs Holt-Winters (same data)

install.packages("wavethresh")
library(wavethresh)
# Original training data
x <- as.numeric(train)

# Pad the series to a power of 2
original_length <- length(x)
padded_length <- 2^ceiling(log2(original_length))

x_padded <- c(
  x,
  rep(tail(x, 1), padded_length - original_length)
)

# Estimate wavelet decomposition
wd_obj <- wavethresh::wd(
  x_padded,
  filter.number = 1,
  family = "DaubExPhase"
)

# Estimate a noise threshold from the finest-scale coefficients
finest_coeff <- wavethresh::accessD(
  wd_obj,
  level = 0
)

sigma_hat <- median(abs(finest_coeff - median(finest_coeff))) / 0.6745

threshold <- sigma_hat * sqrt(2 * log(length(x_padded)))

# Threshold detail coefficients
wd_denoised <- wavethresh::threshold(
  wd_obj,
  policy = "universal",
  type = "soft",
  by.level = FALSE,
  value = threshold
)

# Reconstruct the denoised series
x_denoised <- wavethresh::wr(wd_denoised)

# Remove padded observations
x_denoised <- x_denoised[1:original_length]

# Restore time series structure
train_denoised <- ts(
  x_denoised,
  start = start(train),
  frequency = frequency(train)
)

plot(
  train,
  main = "Original and Wavelet-Denoised Temperature",
  xlab = "Year",
  ylab = "Temperature",
  col = "black"
)

lines(
  train_denoised,
  col = "red",
  lwd = 2
)

legend(
  "topleft",
  legend = c("Original", "Wavelet-denoised"),
  col = c("black", "red"),
  lty = 1,
  bty = "n"
)


# Fit ARIMA to the denoised training data
model_wavelet_arima <- auto.arima(
  train_denoised,
  seasonal = TRUE,
  stepwise = FALSE,
  approximation = FALSE
)

summary(model_wavelet_arima)

# Forecast 12 months
forecast_wavelet_arima <- forecast(
  model_wavelet_arima,
  h = 12
)

library(ggplot2)
autoplot(forecast_wavelet_arima) +
  ggtitle("Wavelet + Box-Jenkins Forecast")

#holt-winter
model_hw <- hw(
  train,
  h = 12,
  seasonal = "additive"
)

autoplot(model_hw) +
  ggtitle("Holt-Winters Forecast")

# Predictions
pred_wavelet <- as.numeric(
  forecast_wavelet_arima$mean
)

pred_hw <- as.numeric(
  model_hw$mean
)

actual <- as.numeric(test)

# Calculate accuracy measures
comparison <- data.frame(
  Method = c(
    "Wavelet + Box-Jenkins",
    "Holt-Winters"
  ),
  
  MAE = c(
    mean(abs(actual - pred_wavelet)),
    mean(abs(actual - pred_hw))
  ),
  
  RMSE = c(
    sqrt(mean((actual - pred_wavelet)^2)),
    sqrt(mean((actual - pred_hw)^2))
  ),
  
  MAPE = c(
    mean(abs((actual - pred_wavelet) / actual)) * 100,
    mean(abs((actual - pred_hw) / actual)) * 100
  )
)

print(comparison)

comparison_plot <- ts(
  cbind(
    Actual = actual,
    Wavelet_ARIMA = pred_wavelet,
    Holt_Winters = pred_hw
  ),
  start = start(test),
  frequency = 12
)

matplot(
  time(comparison_plot),
  comparison_plot,
  type = "l",
  lty = 1,
  col = c("black", "red", "blue"),
  xlab = "Year",
  ylab = "Temperature",
  main = "Actual vs Forecasted Temperature"
)

legend(
  "topleft",
  legend = c(
    "Actual",
    "Wavelet + ARIMA",
    "Holt-Winters"
  ),
  col = c("black", "red", "blue"),
  lty = 1,
  bty = "n"
)


#(a) Wavelet + ARIMA hybrid 
J  <- 3                                      # decomposition levels
wf <- "la8"                                 
mo <- modwt(as.numeric(train), wf = wf, n.levels = J)   

par(mfrow = c(J + 2, 1), mar = c(2, 4, 1, 1))
plot.ts(as.numeric(train), ylab = "Original")
for (nm in names(mo)) plot.ts(mo[[nm]], ylab = nm)
par(mfrow = c(1, 1))

# ARIMA on every component, forecast, then add up
comp_fc <- sapply(names(mo), function(nm) {
  comp <- ts(mo[[nm]], frequency = frequency(train))
  fit  <- auto.arima(comp, seasonal = TRUE)
  as.numeric(forecast(fit, h = h)$mean)
})
wave_fc <- ts(rowSums(comp_fc), start = start(test), frequency = 12)
max(abs(rowSums(sapply(mo, identity)) - as.numeric(train)))


# ---- (b) Holt-Winters -------------------------------------------------------
hw_add <- HoltWinters(train, seasonal = "additive")
hw_mul <- tryCatch(HoltWinters(train, seasonal = "multiplicative"),
                   error = function(e) NULL)
fc_hwa <- forecast(hw_add, h = h)
fc_hwm <- if (!is.null(hw_mul)) forecast(hw_mul, h = h) else NULL
hw_add$alpha; hw_add$beta; hw_add$gamma

# Calculate MAPE while excluding zero actual values
calculate_mape <- function(actual, predicted) {
  actual <- as.numeric(actual)
  predicted <- as.numeric(predicted)
  
  valid <- is.finite(actual) &
    is.finite(predicted) &
    actual != 0
  
  if (!any(valid)) return(NA_real_)
  
  mean(abs((actual[valid] - predicted[valid]) /
             actual[valid])) * 100
}

# ---- (c) Comparison ---------------------------------------------------------
# Calculate accuracy metrics for each model
get_metrics <- function(predicted, actual) {
  predicted <- as.numeric(predicted)
  actual <- as.numeric(actual)
  
  c(
    RMSE = sqrt(mean((actual - predicted)^2, na.rm = TRUE)),
    MAE  = mean(abs(actual - predicted), na.rm = TRUE),
    MAPE = calculate_mape(actual, predicted)
  )
}

res <- rbind(
  "Wavelet+ARIMA"   = get_metrics(wave_fc, test),
  "Holt-Winters(A)" = get_metrics(fc_hwa$mean, test),
  "Seasonal ARIMA"  = get_metrics(fc1$mean, test)
)

if (!is.null(fc_hwm)) {
  res <- rbind(
    res,
    "Holt-Winters(M)" = get_metrics(fc_hwm$mean, test)
  )
}

print(round(res, 4))

cat("MAPE calculated excluding zero actual values.\n")
cat("Number of zero actual values:",
    sum(as.numeric(test) == 0, na.rm = TRUE), "\n")

cat("Best model by RMSE:",
    rownames(res)[which.min(res[, "RMSE"])], "\n")

plot(test, type = "o", ylim = range(test, wave_fc, fc_hwa$mean), main = "Test-set forecasts")
lines(wave_fc,     col = "blue",  type = "o")
lines(fc_hwa$mean, col = "red",   type = "o")
lines(fc1$mean,    col = "darkgreen", type = "o")
legend("topleft", c("Actual", "Wavelet+ARIMA", "Holt-Winters", "SARIMA"),
       col = c("black", "blue", "red", "darkgreen"), lty = 1, pch = 1, cex = 0.8)






#Question-3:
#Data:
rem <- c(1871.49,1810.10,1726.71,1646.87,1553.70,1630.66,1704.53,1494.47,1859.73,
         2010.81,1885.34,1837.27,2096.32,2036.93,1539.60,1525.54,1595.17,1699.70,
         1958.87,1560.48,2022.47,1684.91,1691.66,2199.08,1973.15,1599.45,1334.35,
         1971.43,1930.04,1991.26,2113.15,2164.56,1997.07,2044.23,2254.93,2538.60,
         1913.77,2224.15,2404.11,2395.08,2199.99,2638.78,2185.23,2527.65,3295.63,
         2752.33,2969.56,2822.53,2477.87,2421.89,2685.88)
length(rem) 
#Time series object:
r <- ts(rem, start = c(2021, 7), frequency = 12)
print(r)

library(ggplot2)

remit_df <- data.frame(
  Date = seq(
    from = as.Date("2021-07-01"),
    by = "month",
    length.out = length(rem)
  ),
  Remittance = rem
)

ggplot(remit_df, aes(x = Date, y = Remittance)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Monthly Remittance Inflows to Bangladesh",
    x = "Month",
    y = "Remittance Inflows"
  ) +
  theme_minimal()

#visuailization:
plot(r, main = "Monthly remittance inflow to Bangladesh", ylab = "Million USD",lwd=2,col="blue")
grid()
abline(v = 2024 + 7/12, col = "red", lty = 2)       # August 2024
seasonplot(r, year.labels = TRUE, main = "Seasonal plot")
plot(decompose(r, type = "multiplicative"))

#Variance stabilization
# Seasonal plot
monthplot(
  r,
  main = "Monthly Seasonal Pattern of Remittances",
  ylab = "Remittance"
)

# Autocorrelation
acf(
  r,
  lag.max = 36,
  main = "ACF of Remittance Series"
)

# Partial autocorrelation
pacf(
  r,
  lag.max = 36,
  main = "PACF of Remittance Series"
)

# Classical decomposition
remit_decomp <- decompose(
  r,
  type = "additive"
)

plot(remit_decomp)

# Augmented Dickey-Fuller test
adf.test(r)

# KPSS test
kpss.test(r)


remit_diff <- diff(r)

plot(
  remit_diff,
  main = "First Difference of Remittance Series",
  xlab = "Year",
  ylab = "Change in Remittance"
)

adf.test(remit_diff)
kpss.test(remit_diff)

acf(
  remit_diff,
  main = "ACF of Differenced Remittance"
)

pacf(
  remit_diff,
  main = "PACF of Differenced Remittance"
)

par(mfrow=c(2,2))
acf(r,main="ACF")
pacf(r,main="PACF")
acf(remit_diff,main="ACF diff")
pacf(remit_diff,main="PACF diff")

rl<-log(r)
install.packages("forecast")
library(forecast)
c(
  ARIMA       = forecast::AICc(m_a),
  ARIMA_dummy = forecast::AICc(m_b),
  ETS         = forecast::AICc(m_c)
)

step_dummy <- ts(as.numeric(time(r) >= 2024 + 7/12), start = c(2021, 7), frequency = 12)

m_a <- auto.arima(rl, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
m_b <- auto.arima(rl, xreg = step_dummy, seasonal = TRUE,
                  stepwise = FALSE, approximation = FALSE)
m_c <- ets(rl)
c(ARIMA = AICc(m_a), ARIMA_dummy = AICc(m_b), ETS = AICc(m_c))
summary(m_a); summary(m_b)

# Step 5: hold-out validation (last 6 months) to choose between models
h6 <- 6; nr <- length(r)
tr <- window(rl, end = time(rl)[nr - h6]); te <- window(r, start = time(r)[nr - h6 + 1])
sd_tr <- window(step_dummy, end = time(rl)[nr - h6])
sd_te <- window(step_dummy, start = time(r)[nr - h6 + 1])

f_a <- forecast(auto.arima(tr, stepwise = FALSE), h = h6)
f_b <- forecast(auto.arima(tr, xreg = sd_tr, stepwise = FALSE), h = h6, xreg = sd_te)
f_c <- forecast(ets(tr), h = h6)
cmp <- rbind(ARIMA       = accuracy(exp(f_a$mean), te)[, c("RMSE", "MAE", "MAPE")],
             ARIMA_dummy = accuracy(exp(f_b$mean), te)[, c("RMSE", "MAE", "MAPE")],
             ETS         = accuracy(exp(f_c$mean), te)[, c("RMSE", "MAE", "MAPE")])
print(round(cmp, 3))

best <- m_b; use_x <- TRUE
checkresiduals(best)
shapiro.test(residuals(best))
Box.test(residuals(best), lag = 12, type = "Ljung-Box", fitdf = length(coef(best)))

fut <- ts(rep(1, 6), start = c(2025, 10), frequency = 12)
fc_final <- if (use_x) forecast(best, h = 6, xreg = fut) else forecast(best, h = 6)
fc_tab <- data.frame(
  Month = format(seq(as.Date("2025-10-01"), by = "month", length.out = 6), "%b %Y"),
  Forecast = round(exp(fc_final$mean), 2),
  Lo95 = round(exp(fc_final$lower[, 2]), 2),
  Hi95 = round(exp(fc_final$upper[, 2]), 2))
print(fc_tab)
plot(r, xlim = c(2021.5, 2026.3), ylim = c(1200, 4000), main = "6-month forecast")
lines(exp(fc_final$mean), col = "blue", type = "o")
lines(exp(fc_final$lower[, 2]), col = "grey", lty = 2)
lines(exp(fc_final$upper[, 2]), col = "grey", lty = 2)

# Forecast table (original scale)
fc_tab <- data.frame(
  Month = format(seq(as.Date("2025-10-01"),
                     by = "month", length.out = 6), "%b %Y"),
  Forecast = as.numeric(round(exp(fc_final$mean), 2)),
  Lo95 = as.numeric(round(exp(fc_final$lower[, 2]), 2)),
  Hi95 = as.numeric(round(exp(fc_final$upper[, 2]), 2))
)

print(fc_tab)

# Convert forecast values back to original scale
fc_mean <- exp(fc_final$mean)
fc_lower <- exp(fc_final$lower[, 2])
fc_upper <- exp(fc_final$upper[, 2])

# Time indices
hist_time <- as.numeric(time(r))
fc_time <- as.numeric(time(fc_mean))

# Publication-style plot
par(
  mar = c(5, 5, 4, 2),
  mgp = c(3, 0.8, 0),
  las = 1,
  family = "sans"
)

plot(
  r,
  type = "n",
  xlim = c(2021.5, 2026.3),
  ylim = c(1200, 4000),
  xlab = "Year",
  ylab = "Value",
  main = "Time Series and 6-Month Forecast",
  xaxt = "n",
  bty = "l",
  cex.main = 1.35,
  cex.lab = 1.15,
  cex.axis = 1
)

axis(1, at = 2021:2026, labels = 2021:2026)
abline(h = seq(1200, 4000, by = 400),
       col = "grey90", lty = 1)

# Historical observations
lines(r, col = "black", lwd = 1.8)
points(r, col = "black", pch = 16, cex = 0.65)

# Forecast confidence interval
polygon(
  c(fc_time, rev(fc_time)),
  c(fc_lower, rev(fc_upper)),
  col = adjustcolor("#3B82F6", alpha.f = 0.16),
  border = NA
)

# Confidence interval boundaries
lines(fc_time, fc_lower, col = "grey50", lty = 2, lwd = 1.2)
lines(fc_time, fc_upper, col = "grey50", lty = 2, lwd = 1.2)

# Forecast line and points
lines(fc_time, fc_mean, col = "#1769FF", lwd = 2.5)
points(fc_time, fc_mean, col = "#1769FF", pch = 16, cex = 0.9)

# Historical / forecast divider
abline(v = 2025.75, col = "grey40", lty = 2, lwd = 1.2)

# Annotations
text(2023.3, 3850, "Historical", cex = 1, col = "black")
text(2026.2, 3850, "Forecast (6 months)",
     cex = 0.9, col = "#1769FF", adj = 1)

# Legend
legend(
  "bottom",
  inset = -0.17,
  legend = c("Observed", "Forecast", "95% CI"),
  col = c("black", "#1769FF", adjustcolor("#3B82F6", 0.3)),
  lty = c(1, 1, NA),
  lwd = c(1.8, 2.5, NA),
  pch = c(16, 16, 15),
  pt.cex = c(0.8, 0.9, 1.5),
  horiz = TRUE,
  bty = "n",
  xpd = TRUE,
  cex = 0.95
)






# QUESTION 4 : Stock market data  (last digit = 8)
###############################################################################
# QUESTION 4
install.packages("quantmod")   # only once
library(quantmod)

Q3 <- function(lastDigit){
  sym <- c("JNJ","JPM","XOM","ADBE","META","DIS","V","BRK-B","PG","BAC")
  getSymbols(sym[lastDigit + 1], src = "yahoo",
             from = "2022-01-01", to = "2024-12-31")
}

Q3(8)                    # creates an object called PG in your workspace
PG   <- getSymbols("PG", src = "yahoo", from = "2022-01-01", to = "2024-12-31",
                   auto.assign = FALSE)
data <- PG$PG.Close
str(PG); head(data)

par(mfrow=c(1,1))
# Step 1: explore
px <- ts(as.numeric(data))
plot(px, main = "Closing price", ylab = "Price")
summary(px); sum(is.na(px))

# Step 2: stationarity on price and on log-returns
adf.test(px); kpss.test(px)                    # price: non-stationary expected
ret <- diff(log(px))
plot(ret, main = "Log returns")
adf.test(ret); kpss.test(ret)
tsdisplay(px, main = "Price"); tsdisplay(ret, main = "Log returns")

# Step 3: ARIMA with hold-out (last 10% as test)
N <- length(px); ntest <- round(0.10 * N)
trn <- ts(as.numeric(px)[1:(N - ntest)])
tst <- as.numeric(px)[(N - ntest + 1):N]

fit4 <- auto.arima(trn, stepwise = FALSE, approximation = FALSE)
summary(fit4)
checkresiduals(fit4)
Box.test(residuals(fit4), lag = 20, type = "Ljung-Box", fitdf = length(coef(fit4)))
shapiro.test(residuals(fit4)); jarque.bera.test(residuals(fit4))

# Step 4: volatility check (ARCH effects) and GARCH
ArchTest(residuals(fit4), lags = 12)           # p < 0.05 => ARCH effects
Box.test(residuals(fit4)^2, lag = 12, type = "Ljung-Box")

spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                   distribution.model = "std")
garch_fit <- ugarchfit(spec, data = diff(log(trn)))
show(garch_fit)
plot(garch_fit, which = "all")

# Step 5: forecast + accuracy vs benchmarks
fc4   <- forecast(fit4, h = ntest)
naive <- naive(trn, h = ntest)
drift <- rwf(trn, h = ntest, drift = TRUE)
ets4  <- forecast(ets(trn), h = ntest)
res4 <- rbind(ARIMA = accuracy(fc4,   tst)["Test set", c("RMSE", "MAE", "MAPE")],
              Naive = accuracy(naive, tst)["Test set", c("RMSE", "MAE", "MAPE")],
              Drift = accuracy(drift, tst)["Test set", c("RMSE", "MAE", "MAPE")],
              ETS   = accuracy(ets4,  tst)["Test set", c("RMSE", "MAE", "MAPE")])
print(round(res4, 4))
plot(fc4); lines(ts(tst, start = length(trn) + 1), col = "red")

# Step 6: final model on full data, forecast next 30 observations
final4 <- auto.arima(px, stepwise = FALSE, approximation = FALSE)
plot(forecast(final4, h = 30), main = "30-step-ahead forecast")




