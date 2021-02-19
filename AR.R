library(tidyverse)
library(urca)
library(forecast)
library(tseries) 
library(TSstudio)
gdp <- read_csv("GDP/gdp.csv") %>% ts()

gdp <- gdp1981 %>% ts()

library(forecast)
gdp1 <- gdp[278:468,] %>% select(GDP) %>% ts()

train <- gdp1[1:167,] %>% select(GDP)
test <- gdp1[168:191,] %>% select(GDP)


#####################
h <- 6

split_inf <- ts_split(gdp, sample.out = h)
training <- split_inf$train 
testing <- split_inf$test
length(training)
length(testing)
arima211 <- arima(training, order = c(1,0,0))
autoplot(arima211)
check_res(arima211)

fcast1 <- forecast(arima211, h = h)
test_forecast(actual = gdp, forecast.obj = fcast1, test = testing)

accuracy(fcast1,testing)