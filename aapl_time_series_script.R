# load packages
library(tidyverse) 
library(lubridate) 
library(stats) # predict, acf, pacf, AIC fns
library(moments) # skewness/kurtosis/jarque.test fn
library(zoo) # rollmean fn
library(fUnitRoots) # adfTest fn
library(rugarch) # ugarchspec, ugarchfit fns


# load daily stock price data from yahoo finance
aapl_daily <- read.csv("~/Econometrics/AAPL.csv")

################################################ 
# SCRUB
################################################

# convert daily data into monthly data 
aapl_monthly <- aapl_daily %>%
  mutate(Date = ymd(Date), # convert Date to date format
         year = year(Date), # create year and month vectors
         month = month(Date)) %>% 
  group_by(year, month) %>% # group by year and month
  arrange(Date) %>% # arrange/sort by month desc
  filter(row_number() == max(row_number())) # select last row of month 


aapl_monthly[,8:9] = NULL # drop year and month vectors that were just created above

# difference monthly data
aapl_monthly$dclose <- c(NA, 100*diff(log(aapl_monthly$Adj.Close))) # create dclose variable, difference log variable to calculate returns
aapl_monthly <- aapl_monthly[-1,] # drop NA row due to difference


#############################################
# EXPLORE
#############################################

# summary statistics
aapl_monthly_stats <- aapl_monthly %>%
  select(dclose) %>%
  summarise(
    mean = mean(dclose), 
    median = median(dclose), 
    max = max(dclose), 
    min = min(dclose), 
    sd = sd(dclose), 
    skewness = skewness(dclose), # 3rd moment
    kurtosis = kurtosis(dclose), # 4th moment
    count = length(dclose)
  )

aapl_monthly_stats # call monthly stats


# create histogram
ggplot(aapl_monthly, aes(x = dclose)) +
  geom_histogram(breaks = seq(min(aapl_monthly$dclose), max(aapl_monthly$dclose)+10, 5), 
                 color = "black") 

# formal test for normality
jarque.test(aapl_monthly$dclose)


# create criterion tables
aic_table <- array(NA, c(4,4,2)) # (ar, ma, tables | index 1)
for (ar in 0:3) { # index 0, should be 1 less than above
  for (ma in 0:3) {
    arma <- arima(aapl_monthly$dclose[aapl_monthly$Date <= "2019-12-31"], order = c(ar, 0, ma)) # estimate model iteratively
    aic_table[ar+1, ma+1, 1] <- AIC(arma) # AIC 
    aic_table[ar+1, ma+1, 2] <- AIC(arma, k = log(nrow(aapl_monthly))) # SBIC
  }
}

aic_table
which.min(aic_table[,,1]) # AIC 
which.min(aic_table[,,2]) # SBIC 


# check correlograms for autocorrelation process
par(mfrow=c(1,2))
acf(aapl_monthly$dclose, lag.max = 36, main = "Autocorrelation Function")
pacf(aapl_monthly$dclose, lag.max = 36, main = "Partial Autocorrelation Function")

# test for unit root
adfTest(aapl_monthly$dclose, lags = 36, type = "c")

adfgls = urersTest(aapl_monthly$dclose, type = "DF-GLS", model = "trend", lag.max = 36, doplot = FALSE)
adfgls@test$test@teststat
adfgls@test$test@cval

# formal test for white noise process
tsdiag(arima(aapl_monthly$dclose, order = c(0, 0, 0)), gof.lag = 36) 
Box.test(aapl_monthly$dclose, lag = 36, fitdf = 0, type = "Lj") # ljung-box test



################################################ 
# MODEL
################################################

# fit/estimate arma00 model
arma00 <- arima(aapl_monthly$dclose[aapl_monthly$Date <= "2019-12-31"], order = c(0, 0, 0)) # fit arma00 model

# fixed sample size/dynamic sample 1 step ahead forecast
arma00_fc_fixed <- predict(arma00, n.ahead = 19)$pred # dynamic/fixed sample
arma00_fc_fixed <- as_tibble(arma00_fc_fixed) # turn estimates into tibble

arma00_fc_fixed <- arma00_fc_fixed %>%
  rename(arma00_fc = x) %>% # rename the predicted values
  mutate(Date = as.Date(aapl_monthly$Date[aapl_monthly$Date > "2019-12-31"]), # bring in out of sample date field from aapl_monthly
         residual = as.numeric(aapl_monthly$dclose[aapl_monthly$Date > "2019-12-31"] - arma00_fc), # calculate residuals
         residual2 = as.numeric(residual^2) # calculate residuals squared
  ) %>% 
  select(Date, arma00_fc, residual, residual2)  # reorder fields to match data layout

# static sample size/rolling sample 1 step ahead forecast
arma00_fc_rolling <- aapl_monthly %>%
  mutate(fc = as.numeric(c(rep(NA, 372), rollmean(aapl_monthly$dclose, 373, align = "right"))), # create vector of NAs and predicted arma00 (mean)
         residual = as.numeric(dclose - fc), # calculate the residual 
         residual2 = as.numeric(residual^2)) %>% # calculate the residual squared
  filter(!is.na(fc)) %>% # filter out the NAs
  mutate(arma00_fc = as.ts(fc)) %>%
  select(Date, arma00_fc, residual, residual2) # select Date, forecasted values and residuals


# fit/estimate arbitrary arma11 model for comparison
arma11 <- arima(aapl_monthly$dclose[aapl_monthly$Date <= "2019-12-31"], order = c(1, 0, 1)) # fit arma11

arma11_fc_fixed <- predict(arma11, n.ahead = 19)$pred # dynamic / fixed sample

# fixed sample size/dynamic sample 1 step ahead forecast
arma11_fc_fixed <- as_tibble(arma11_fc_fixed) %>% # turn estimates into tibble
  rename(arma11_fc = x) %>% # rename predicted values
  mutate(Date = as.Date(aapl_monthly$Date[aapl_monthly$Date > "2019-12-31"]), # bring in out of sample date field from aapl_monthly
         residual = as.numeric(aapl_monthly$dclose[aapl_monthly$Date > "2019-12-31"] - arma11_fc), # calculate the residuals
         residual2 = as.numeric(residual^2) # calculate the residuals squared
  ) %>%
  select(Date, arma11_fc, residual, residual2) # reorder fields to match data layout


#####################################
# EVALUATE
#####################################

# calculate mean squared errors (mse)
arma00_fc_fixed_mse <- mean(arma00_fc_fixed$residual2) # calculate the fixed sample arma00 model's out of sample mean squared error
arma00_fc_rolling_mse <-  mean(arma00_fc_rolling$residual2) # calculate the rolling sample arma00 model's out of sample mean squared error
arma11_fc_fixed_mse <- mean(arma11_fc_fixed$residual2) # calculate the fixed sample arma11 model's out of sample mean squared error


mse <- matrix(c("arma00_fc_fixed", 
                arma00_fc_fixed_mse,
                "arma00_fc_rolling", 
                arma00_fc_rolling_mse,
                "arma11_fc_fixed", 
                arma11_fc_fixed_mse),
              2, 3)
mse

# view predicted returns compared to actuals 
par(mfrow=c(1,1))
plot(aapl_monthly$Date[372:391], aapl_monthly$dclose[372:391], type = "l", ylab = "Percentage Returns", xlab = "Time")
lines(arma11_fc_fixed$Date, arma11_fc_fixed$arma11_fc, lwd = 7, col = "green")
lines(arma00_fc_rolling$Date, arma00_fc_rolling$arma00_fc, lwd = 4, col = "blue")
lines(arma00_fc_fixed$Date, arma00_fc_fixed$arma00_fc, lwd = 3, col = "red")
legend("topright", legend = c("arma00_fc_fixed", "arma00_fc_rolling", "arma11_fc_fixed"), lwd = 3, col = c("red", "blue", "green"))


