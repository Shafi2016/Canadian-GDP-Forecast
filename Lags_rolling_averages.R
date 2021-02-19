library(tidymodels)
library(modeltime)
library(prophet)
library(tidyquant)
# Core 
library(tidyverse)
library(lubridate)
library(timetk)
library(plotly)




gt2 <- data %>%
  tk_augment_lags(.value = c(ALTSALES : PAYEMS), .lags = c(1:5)) %>%
  
  tk_augment_slidify(
    .value  = c(GDP :HOUST),
    .period  = c(2:10),
    .f       = mean,
    .align   = "center", 
    .partial = TRUE)    %>% 
  
    tk_augment_slidify(
       .value  = c(GDP :IP_new ),
       .period  = c(2:10),
       .f       = mean,
       .align   = "center",
       .partial = TRUE) %>% drop_na()


write.csv(gt2, file = 'canus.csv')