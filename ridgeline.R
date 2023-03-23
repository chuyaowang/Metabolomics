library(ggplot2)
library(openxlsx)
library(dplyr)
library(future.apply)

data <- read.xlsx("data/ninefigures.xlsx")
time <- round(data$time,digits = 3)
ints <- round(data[,2:ncol(data)])

df <- future_lapply(ints, function(x){ # x: each column in ints
  temp <- c()
  for (i in seq_along(x)) {
    temp <- append(temp,rep(time[i],times=x[i]))
  }
  return(temp)
})