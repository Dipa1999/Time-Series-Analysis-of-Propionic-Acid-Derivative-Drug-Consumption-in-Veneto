# Consumption of propionic and derivative drugs in the Veneto region of Italy

# Set the working directory to the same path as the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Load libraries
library(tidyverse)
library(dplyr)
library(conflicted)
library(forecast)
library(lmtest) 
library(prophet)

conflict_prefer('filter', 'dplyr')


# Data available at: https://www.aifa.gov.it/spesa-e-consumo-relativi-al-flusso-della-farmaceutica-convenzionata-e-degli-acquisti-diretti

# Load data

data <- bind_rows(read.csv('dati2016_100519.csv', 
                           sep = '|'), 
                  read.csv('dati2017_100519.csv', 
                           sep = '|'), 
                  read.csv('dati2018_23.09.2020.csv', 
                           sep = '|'), 
                  read.csv('dati2019_23.09.2020.csv', 
                           sep = '|'), 
                  read.csv('dati2020_22.10.2021.csv', 
                           sep = '|'), 
                  read.csv('dati2021_24.10.2022.csv', 
                           sep = '|'), 
                  read.csv('dati2022_07.02.2024.csv', 
                           sep = '|'), 
                  read.csv('dati2023_15.01.2025.csv', 
                           sep = '|'))




# Select propionic acid-based drugs and charged to the caregiver

library(dplyr)

data <- data %>% 
  filter(regione == 'VENETO', 
         atc4 == 'M01AE') %>%
  mutate(n_packs = numero_confezioni_traccia + numero_confezioni_convenzionata) %>%  # All packs
  select(anno, mese, n_packs) %>%  
  group_by(anno, mese) %>%  
  summarise(n_packs = sum(n_packs, na.rm = TRUE)) %>%
  mutate(data = make_date(anno, mese, 1)) %>% 
  rename(year = anno, 
         month = mese)


h <- 12  # Monthly frequency
ts_data <- ts(data$n_packs, 
              start = c(min(data$year), min(data$month)), 
              frequency = h)  # Monthly frequency


plot(ts_data, 
     type = 'b', 
     xlab = 'Year', 
     ylab = 'Number of packs', 
     main = 'Consumption of propionic acid-based drugs in the Veneto region of Italy')  # Remove x-axis labels
# Add a grid
grid()


# Plot ACF and PACF
par(mfrow = c(1, 2))
acf(ts_data, 
    lag.max = 24, 
    main = '')
pacf(ts_data,
     lag.max = 24, 
     main = '')
par(mfrow = c(1, 1))


# Seasonal plot
seasonplot(ts_data, 
           ylab="Number of packs", 
           xlab="Months", 
           main="Seasonal plot", 
           year.labels=T, 
           year.labels.left=T, 
           col=1:8, pch=19)
# Add a grid
grid()


# We can see that the series has a trend and a seasonality. 


# EDA
# Check for trend
t1 <- 1:nrow(data)
lm1 <- lm(data$n_packs ~ t1)
summary(lm1)


ggplot(data, 
       aes(x = data, 
           y = n_packs)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = 'lm',
              col='red') +  # Add trend line
  scale_x_date(date_breaks = '1 year', 
               date_labels = '%Y') +  # Change x-axis labels
  labs(x = 'Year', 
       y = 'Number of packs', 
       title = 'Trend') +
  theme_minimal()

# We can see the covid and lockdowns effect on the series


# Analyse where the are the strokes

# Evaluate the distribution of the residuals on trend
# evaluating until 2020

ts_0 <- window(ts_data, end = c(2019, 12))
time_0 <- 1:length(ts_0)
fit_0 <- lm(ts_0 ~ time_0)

sd_0 <- sd(residuals(fit_0))
alpha <- 0.05
low_99 <- qt(alpha / 2, df = length(ts_0) - 2) * sd_0 
up_99 <- qt(1 - alpha / 2, df = length(ts_0) - 2) * sd_0

l_01 <- ts(fit_0$coefficients[1] + 
             (fit_0$coefficients[2] * 1:length(ts_data)) + 
             low_99, 
           start = c(2016, 1), 
           frequency = h)
u_99 <- ts(fit_0$coefficients[1] + 
             (fit_0$coefficients[2] * 1:length(ts_data)) + 
             up_99, 
           start = c(2016, 1), 
           frequency = h)

plot(ts_data, 
     ylim = c(min(c(ts_data, l_01)), 
              max(c(ts_data, u_99))),
     type = 'b', 
     xlab = 'Year',
     ylab = 'Number of packs', 
     main = 't-based Outliers Detection')
lines(l_01, col = 'red')
lines(u_99, col = 'blue')
grid()
legend('topright', 
       legend = c('Observed data', 
                  '95% CI lower bound', 
                  '95% CI upper bound'), 
       col = c('black', 
               'red', 
               'blue'), 
       lty = 1, 
       cex = 0.8)


# Outliers date
as.yearmon(time(ts_data)[which(ts_data < l_01 | ts_data > u_99)])

# We can set January 2020 as outliers and manage April and May 2020 as lockdowns



### Introducing lockdowns dummy
# covid_ts <- ifelse(data$year == 2020 & 
#                      data$month == 4,   # One month later lockdown
#                    1, 0) %>% 
#   ts(start = c(min(data$year), min(data$month)),
#      frequency = h)

# Covid awareness
covid_awareness <- ifelse(data$year == 2020 & 
                            data$month == 1,   # Jan 2020
                          1, 0) %>% 
  ts(start = c(min(data$year), min(data$month)),
     frequency = h)

lockdown_ts <- ifelse(data$year == 2020 & 
                        data$month %in% c(4, 5), 
                      1, 0) %>% 
  ts(start = c(min(data$year), min(data$month)),
     frequency = h)


# Split the ts in train and test, train until 2022
ts_train <- window(ts_data, end = c(2022, 12))
# covid_train <- window(covid_ts, end = c(2022, 12))
awareness_train <- window(covid_awareness, end = c(2022, 12))
lockdown_train <- window(lockdown_ts, end = c(2022, 12))

ts_test <- window(ts_data, start = c(2023, 1))
# covid_test <- window(covid_ts, start = c(2023, 1))
awareness_test <- window(covid_awareness, start = c(2023, 1))
lockdown_test <- window(lockdown_ts, start = c(2023, 1))




# Fit Models --------------------------------------------------------------



## ARIMA
# without covid variables
fit_arima <- auto.arima(ts_train)
summary(fit_arima)

# Check residuals
checkresiduals(fit_arima)
acf(residuals(fit_arima))
pacf(residuals(fit_arima))
Box.test(residuals(fit_arima), lag = 24, type = 'Ljung-Box')

# No significant autocorrelation in the residuals

# Setting the layout for the plots
layout(matrix(c(1, 1,
                2, 3), 
              nrow = 2, byrow = TRUE))
plot(resid(fit_arima), 
     main = 'Residuals of ARIMA(0,1,1)(1,0,0)[12] ', 
     ylab = 'Residuals', 
     xlab = 'Time')
acf(residuals(fit_arima), 
    lag.max = 24, 
    ylab = 'ACF', 
    main = '')
pacf(residuals(fit_arima),
     lag.max = 24, 
     ylab = 'PACF', 
     main = '')
layout(1)  # Reset layout


# Forecast
forecast_arima <- forecast(fit_arima, h = h)
plot(forecast_arima)
lines(ts_test, col = 'red')



## ARIMA with covid's lockdowns dummy
fit_arima_covid <- auto.arima(ts_train, 
                              xreg = cbind(
                                # covid = covid_train, 
                                awareness = awareness_train,
                                lock = lockdown_train))
summary(fit_arima_covid)

# Check residuals
checkresiduals(fit_arima_covid)
acf(residuals(fit_arima_covid))
pacf(residuals(fit_arima_covid))
Box.test(residuals(fit_arima_covid), lag = 24, type = 'Ljung-Box')

# No significant autocorrelation in the residuals

# Setting the layout for the plots
layout(matrix(c(1, 1,
                2, 3), 
              nrow = 2, byrow = TRUE))
plot(resid(fit_arima_covid),
     main = 'Residuals of ARIMA(0,1,1)(1,0,1)[12] ', 
     ylab = 'Residuals', 
     xlab = 'Time')
acf(residuals(fit_arima_covid),
    lag.max = 24, 
    ylab = 'ACF', 
    main = '')
pacf(residuals(fit_arima_covid),
     lag.max = 24, 
     ylab = 'PACF', 
     main = '')
layout(1)  # Reset layout

# Forecast
forecast_arima_covid <- forecast(fit_arima_covid, 
                                 xreg = cbind(
                                   # covid = covid_test, 
                                   awareness = awareness_test, 
                                   lock = lockdown_test))
plot(forecast_arima_covid)
lines(ts_test, col = 'red')
lines(forecast_arima$mean, col = 'green')


# Plot different models
plot(forecast_arima$mean, 
     ylim = c(min(c(forecast_arima$mean, forecast_arima_covid$mean)), 
              max(c(forecast_arima$mean, forecast_arima_covid$mean))))
lines(forecast_arima_covid$mean, col = 'red')



## ARIMAX
season_train <- seasonaldummy(ts_train)

fit_arimax <- auto.arima(ts_train, 
                         xreg = cbind(
                           time = 1:length(ts_train),
                           seasonaldummy(ts_train),
                           # covid = covid_train,
                           awareness = awareness_train,
                           lock = lockdown_train) %>% 
                           `colnames<-`(sub('seasonaldummy\\(ts_train)\\.', '', 
                                            colnames(.)))) # Keep only the month
summary(fit_arimax) # ARIMA(0, 0, 0), need to fix

# Check residuals
checkresiduals(fit_arimax)
acf(residuals(fit_arimax))
pacf(residuals(fit_arimax))
Box.test(residuals(fit_arimax), lag = 24, type = 'Ljung-Box')
# Significant autocorrelation in the residuals

library(tseries)
adf.test(residuals(fit_arimax))
# Significant of non stationary residuals



fit_arimax <- auto.arima(ts_train, 
                         xreg = cbind(
                           time = 1:length(ts_train),
                           seasonaldummy(ts_train), 
                           # covid = covid_train,
                           awareness = awareness_train,
                           lock = lockdown_train) %>% 
                           `colnames<-`(sub('seasonaldummy\\(ts_train)\\.', '', 
                                            colnames(.))), 
                         stationary = F,
                         # max.order = 12, 
                         stepwise = F, 
                         trace = T)
summary(fit_arimax)

# Check residuals
checkresiduals(fit_arimax)
acf(residuals(fit_arimax))
pacf(residuals(fit_arimax))
Box.test(residuals(fit_arimax), lag = 24,  type = 'Ljung-Box')
adf.test(residuals(fit_arimax))
# Significant autocorrelation in the residuals

# We need to add a difference term to make the residuals stationary

# To check the order of difference
data_lm <- cbind(n_packs = ts_data, 
                 time = c(1:length(ts_data)), 
                 seasonaldummy(ts_data), 
                 # covid = covid_ts,
                 awareness = covid_awareness,
                 lockdown = lockdown_ts) %>% 
  as.data.frame() %>% 
  `colnames<-`(sub('seasonaldummy\\(ts_data)\\.', '', 
                   colnames(.)))

train_lm <- data_lm[1:length(ts_train), ]
test_lm <- data_lm[(length(ts_train) + 1) : nrow(data_lm), ]

fit_lm <- lm(n_packs ~ ., 
             data = train_lm)
summary(fit_lm)

resid_lm <- residuals(fit_lm)

# Create ARIMA model on residuals
fit_arima_lm <- Arima(resid_lm, 
                      order = c(0, 1, 1), 
                      seasonal = list(order = c(0, 0, 0), 
                                      period = h))
summary(fit_arima_lm)
checkresiduals(fit_arima_lm)
acf(resid(fit_arima_lm))
pacf(resid(fit_arima_lm))
Box.test(resid(fit_arima_lm), lag = 24, type = 'Ljung-Box')
adf.test(resid(fit_arima_lm))
# Great

# Setting the layout for the plots
layout(matrix(c(1, 1,
                2, 3), 
              nrow = 2, byrow = TRUE))
plot(ts(resid(fit_arima_lm), 
        start = c(min(data$year), min(data$month)), 
        frequency = h),
     main = 'Residuals of ARIMAX(0,1,1)(0,0,0)[12] ', 
     ylab = 'Residuals', 
     xlab = 'Time')
acf(residuals(fit_arima_lm),
    lag.max = 24, 
    ylab = 'ACF', 
    main = '')
pacf(residuals(fit_arima_lm),
     lag.max = 24, 
     ylab = 'PACF', 
     main = '')
layout(1)  # Reset layout


# Fitted values

fitted_arimax <- fitted(fit_lm) + 
  fitted(fit_arima_lm)
fitted_arimax <- ts(fitted_arimax, 
                    start = c(2016, 1), 
                    frequency = h)

# Forecast

forecast_arimax <- predict(fit_lm, 
                           newdata = test_lm) + 
  forecast(fit_arima_lm, 
           h = h)$mean
forecast_arimax <- ts(forecast_arimax, 
                      start = c(2023, 1), 
                      frequency = h)


# Plot fitted values on the real data
plot(ts_train, type = 'b', 
     xlab = 'Year',
     ylab = 'Number of packs', 
     main = 'Fitted values')
lines(fitted(fit_arima), col = 'blue')
lines(fitted(fit_arima_covid), col = 'purple')
lines(fitted_arimax, col = 'green')
legend('topright', 
       legend = c('Real data', 
                  'ARIMA', 
                  'ARIMA with COVID', 
                  'ARIMAX'), 
       col = c('black', 
               'blue', 
               'purple', 
               'green'), 
       lty = 1, 
       cex = 0.8)





## Prophet

df <- data.frame(
  ds = as.Date(ts_train),
  y = as.numeric(ts_train)
)


# Create a dataframe with the holidays

# covid <- data.frame(
#   holiday = 'COVID',
#   ds = as.Date('2020-04-01'),  # COVID
#   lower_window = 0,
#   upper_window = 0,
#   impact = 1  # Positive impact
# )

awareness <- data.frame(
  holiday = 'Awareness',
  ds = as.Date('2020-01-01'),  # awareness
  lower_window = 0,
  upper_window = 0
)

lockdown <- data.frame(
  holiday = 'Lockdown',
  ds = as.Date(c('2020-04-01', '2020-05-01')),  # lockdown
  lower_window = 0,
  upper_window = 0
)


# Create the Prophet model including the holidays
fit_prophet <- prophet(df, 
                       yearly.seasonality = T, 
                       daily.seasonality = F,
                       weekly.seasonality = F,
                       holidays = rbind(
                         # covid, 
                         awareness,
                         lockdown))

summary(fit_prophet)

# Forecast
future <- make_future_dataframe(fit_prophet, 
                                periods = h, 
                                freq = 'month', 
                                include_history = T)
forecast_prophet <- predict(fit_prophet, future)

tail(forecast_prophet[c('ds', 'yhat', 'yhat_lower', 'yhat_upper')])


# Parameters
prophet_plot_components(fit_prophet, forecast_prophet)

# Prediction plot

plot(fit_prophet, forecast_prophet, 
     xlab = 'Year',
     ylab = 'Number of packs')


# Setting the layout for the plots
resid_prophet <- ts_train - predict(fit_prophet)$yhat
layout(matrix(c(1, 1,
                2, 3), 
              nrow = 2, byrow = TRUE))
plot(resid_prophet,
     main = 'Residuals of prophet ', 
     ylab = 'Residuals', 
     xlab = 'Time')
acf(resid_prophet,
    lag.max = 24, 
    ylab = 'ACF', 
    main = '')
pacf(resid_prophet,
     lag.max = 24, 
     ylab = 'PACF', 
     main = '')
layout(1)  # Reset layout


# Dynamic plot
dyplot.prophet(fit_prophet, forecast_prophet, 
               main = 'Forecasted values',
               xlab = 'Year',
               ylab = 'Number of packs') 

# Plot with change points

plot(fit_prophet, forecast_prophet) + 
  add_changepoints_to_plot(fit_prophet, threshold=0)

# Dates corresponding to change points
fit_prophet$changepoints




# Compare models ----------------------------------------------------------





# Plot fitted values on the real data
ts_prophet <- ts(forecast_prophet$yhat, 
                 start = c(min(data$year), min(data$month)), 
                 frequency = h)

ts_prophet_2022 <- window(ts_prophet, end = c(2022, 12))
ts_prophet_2023 <- window(ts_prophet, start = c(2023, 01))


# Plot fitted values on the real data
plot(ts_train, type = 'b', 
     xlab = 'Year',
     ylab = 'Number of packs', 
     main = 'Fitted values')
grid()
lines(fitted(fit_arima), col = 'blue')
lines(fitted(fit_arima_covid), col = 'purple')
lines(fitted_arimax, col = 'green')
lines(ts_prophet_2022, col = 'red')
legend('bottomlef', 
       legend = c('Real data', 
                  'SARIMA', 
                  'SARIMA with COVID', 
                  'ARIMAX', 
                  'Prophet'), 
       col = c('black', 
               'blue', 
               'purple', 
               'green', 
               'red'), 
       lty = 1, 
       cex = 0.8)


# Plot forecasted values
plot(ts_test, type = 'b', 
     xaxt = 'n', 
     ylim = c(min(c(ts_test, 
                    forecast_arima$mean, 
                    forecast_arima_covid$mean, 
                    forecast_arimax, 
                    ts_prophet_2023)), 
              max(c(ts_test, 
                    forecast_arima$mean, 
                    forecast_arima_covid$mean, 
                    forecast_arimax, 
                    ts_prophet_2023)) + 
                2000), 
     xlab = 'Year',
     ylab = 'Number of packs',
     main = 'Forecasted values 2023')
axis(1, at = time(ts_test), labels = month.abb)
grid()
lines(forecast_arima$mean, col = 'blue')
lines(forecast_arima_covid$mean, col = 'purple')
lines(forecast_arimax, col = 'green')
lines(ts_prophet_2023, col = 'red')
legend('bottomlef', 
       legend = c('Real data', 
                  'SARIMA', 
                  'SARIMA with COVID', 
                  'ARIMAX', 
                  'Prophet'), 
       col = c('black', 
               'blue', 
               'purple', 
               'green', 
               'red'), 
       lty = 1, 
       cex = 0.8)





# Validate models ---------------------------------------------------------


## Train

# Create list of residuals' models
residuals_train <- list(resid(fit_arima), 
                        resid(fit_arima_covid), 
                        fitted_arimax - ts_train, 
                        ts_prophet_2022 - ts_train)

# ME (Mean Error)
ME_train <- sapply(residuals_train, mean)

# MAE (Mean Absolute Error)
MAE_train <- sapply(residuals_train, 
                    function(x) mean(abs(x)))

# RMSE (Root Mean Squared Error)
RMSE_train <- sapply(residuals_train, 
                     function(x) sqrt(mean(x^2)))

# MPE (Mean Percentage Error)
MPE_train <- sapply(residuals_train, 
                    function(x) mean((x / ts_train) * 100))

# MAPE (Mean Absolute Percentage Error)
MAPE_train <- sapply(residuals_train, 
                     function(x) mean(abs(x / ts_train) * 100))


## Test

# Create list of residuals' models
# We don't use AIC or BIC since is not available for prophet
residuals_test <- list(forecast_arima$mean - ts_test, 
                       forecast_arima_covid$mean - ts_test, 
                       forecast_arimax - ts_test, 
                       ts_prophet_2023 - ts_test)

# ME (Mean Error)
ME_test <- sapply(residuals_test, mean)

# MAE (Mean Absolute Error)
MAE_test <- sapply(residuals_test, 
                   function(x) mean(abs(x)))

# RMSE (Root Mean Squared Error)
RMSE_test <- sapply(residuals_test, 
                    function(x) sqrt(mean(x^2)))

# MPE (Mean Percentage Error)
MPE_test <- sapply(residuals_test, 
                   function(x) mean((x / ts_test) * 100))

# MAPE (Mean Absolute Percentage Error)
MAPE_test <- sapply(residuals_test, 
                    function(x) mean(abs(x / ts_test) * 100))

evaluation_train <- data.frame(
  models = c('SARIMA', 
             'SARIMA with COVID', 
             'ARIMAX', 
             'Prophet'), 
  ME = round(ME_train, 2),
  MAE = round(MAE_train, 2),
  RMSE = round(RMSE_train, 2),
  MPE = round(MPE_train, 2),
  MAPE = round(MAPE_train, 2)
)


evaluation_test <- data.frame(
  models = c('SARIMA', 
             'SARIMA with COVID', 
             'ARIMAX', 
             'Prophet'), 
  ME = round(ME_test, 2),
  MAE = round(MAE_test, 2),
  RMSE = round(RMSE_test, 2),
  MPE = round(MPE_test, 2),
  MAPE = round(MAPE_test, 2)
)


evaluation_train
evaluation_test

# For latex
library(xtable)
xtable(evaluation_train)
xtable(evaluation_test)

# All models works well, but the prophet model is the best one on train
# and ARIMAX on test

# Calculate the surplus in 2023
surplus <- data.frame(
  models = c('SARIMA', 
             'SARIMA with COVID', 
             'ARIMAX', 
             'Prophet'), 
  surplus = rbind(forecast_arima$mean - ts_test, 
                  forecast_arima_covid$mean - ts_test, 
                  forecast_arimax - ts_test, 
                  ts_prophet_2023 - ts_test)
)

colnames(surplus) <- c('models', 
                       'Jan', 'Feb', 'Mar', 'Apr', 
                       'May', 'Jun', 'Jul', 'Aug', 
                       'Sep', 'Oct', 'Nov', 'Dec')

surplus <- rowsum(surplus[, -1], 
                  group = surplus$models)
rowSums(surplus)
cbind(surplus, rowSums(surplus)) # Total surplus for each model

sum(ts_test)

rowSums(surplus) / sum(ts_test) * 100

surplus_divided <- apply(surplus, 1, 
                         function(row) row / ts_test * 100)
cbind(surplus_divided, rowSums(surplus_divided))
rbind(surplus_divided, colSums(surplus_divided))


# We can see that ARIMA models are more conservative and 
# with less variance than others,
# while ARIMAX and Prophet models perform better in terms of surplus





# Retrain models on the whole dataset -------------------------------------



## ARIMA
final_arima <- Arima(ts_data, 
                     order = c(0, 1, 1), 
                     seasonal = list(order = c(1, 0, 0), 
                                     period = h))
summary(final_arima)

# check residuals
checkresiduals(final_arima)
acf(residuals(final_arima))
pacf(residuals(final_arima))
Box.test(residuals(final_arima), lag = 24, type = 'Ljung-Box')
adf.test(residuals(final_arima))

# Forecast
f_arima <- forecast(final_arima, h = h)



## ARIMA with covid's lockdowns dummy
final_arima_covid <- Arima(ts_data, 
                           order = c(0, 1, 1), 
                           seasonal = list(order = c(1, 0, 1), 
                                           period = h), 
                           xreg = cbind(
                             # covid = covid_ts, 
                             awareness = covid_awareness,
                             lock = lockdown_ts))
summary(final_arima_covid)

# Check residuals
checkresiduals(final_arima_covid)
acf(residuals(final_arima_covid))
pacf(residuals(final_arima_covid))
Box.test(residuals(final_arima_covid), lag = 24, type = 'Ljung-Box')


# Forecast
cov_test <- aw_test <- lock_test <- ts(rep(0, 12), 
                                       start = c(2024, 1), 
                                       frequency = h)
f_arima_covid <- forecast(final_arima_covid, 
                          xreg = cbind(
                            # covid = cov_test, 
                            awareness = aw_test,
                            lock = lock_test))



## ARIMAX
test2_lm <- test_lm
test2_lm$time <- test2_lm$time + 12

final_lm <- lm(n_packs ~ ., 
               data = data_lm)
summary(final_lm)

resid_final <- residuals(final_lm)


# Create ARIMA model on residuals
final_arima_lm <- Arima(resid_final, 
                        order = c(0, 1, 1), 
                        seasonal = list(order = c(0, 0, 0), 
                                        period = h))
summary(final_arima_lm)
checkresiduals(final_arima_lm)
acf(resid(final_arima_lm))
pacf(resid(final_arima_lm))
Box.test(resid(final_arima_lm), lag = 24, type = 'Ljung-Box')
# Little autocorrelation in the residuals
adf.test(resid(final_arima_lm))
# Great



# Forecast
# We keep the same order of the variables
f_arimax <- predict(final_lm, 
                    newdata = test2_lm, 
                    interval = 'prediction', # for estimate each individual month and not the average
                    level = 0.95) + 
  matrix(rep(forecast(final_arima_lm, 
                      h = h)$mean, 3),
         ncol = 3, 
         byrow = T)
f_arimax <- ts(f_arimax, 
               start = c(2024, 1), 
               frequency = h)




## Prophet

df <- data.frame(
  ds = as.Date(ts_data),
  y = as.numeric(ts_data)
)



# Create the Prophet model including the holidays
final_prophet <- prophet(df, 
                         yearly.seasonality = T, 
                         daily.seasonality = F,
                         weekly.seasonality = F,
                         holidays = rbind(
                           # covid, 
                           awareness,
                           lockdown))

# summary(final_prophet)

# Forecast
future <- make_future_dataframe(final_prophet, 
                                periods = h, 
                                freq = 'month', 
                                include_history = T)
f_prophet <- predict(final_prophet, future)

tail(f_prophet[c('ds', 'yhat', 'yhat_lower', 'yhat_upper')])

#prediction plot

plot(final_prophet, f_prophet)

ts_prophet_final <- ts(f_prophet$yhat, 
                       start = c(min(data$year), min(data$month)), 
                       frequency = h)

ts_prophet_2023_b <- window(ts_prophet_final, end = c(2023, 12))
ts_prophet_2024 <- window(ts_prophet_final, start = c(2024, 01))



### Plot forecasted values

# Date start plot
d_start = 2021


plot(window(ts_data, 
            start = d_start), 
     type = 'b', 
     ylim = c(min(c(ts_data[(5*12) : length(ts_data)], f_arima$mean, 
                    f_arima_covid$mean, f_arimax, 
                    f_prophet$yhat[(5*12) : length(f_prophet$yhat)])), 
              max(c(ts_data[(5*12) : length(ts_data)], f_arima$mean, 
                    f_arima_covid$mean, f_arimax, 
                    f_prophet$yhat[(5*12) : length(f_prophet$yhat)]))),  
     xaxt = 'n',
     xlim = c(d_start, 2025), 
     xlab = 'Year',
     ylab = 'Number of packs',
     main = 'Forecasted values 2024')
grid()
# Add years to x-axis
axis(1, 
     at = seq(d_start, 2025, 1), 
     labels = seq(d_start, 2025, 1))
lines(f_arima$mean, col = 'blue')
lines(f_arima_covid$mean, col = 'purple')
lines(f_arimax[, 1], col = 'green')   # mean
lines(f_arimax[, 2], 
      lty = 'dashed', 
      col = 'green')   # lower 0.95
lines(f_arimax[, 3], 
      lty = 'dashed', 
      col = 'green')   # upper 0.95
lines(ts_prophet_2024, col = 'red')
legend('bottomleft', 
       legend = c('Real data', 
                  'SARIMA', 
                  'SARIMA with COVID', 
                  'ARIMAX', 
                  '95% CI ARIMAX',
                  'Prophet'), 
       col = c('black', 
               'blue', 
               'purple', 
               'green', 
               'green',
               'red'), 
       lty = c(1, 
               1, 
               1, 
               1, 
               2,
               1),
       cex = 0.8)


#### END





