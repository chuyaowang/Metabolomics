library(dplyr)
library(ggplot2)
library(readr)

data <- read.csv("rfApp/data/train.csv")
params <- data %>% 
  group_by(y) %>% 
  summarize(mean_x1=mean(x1),sd_x1=sd(x1),mean_x2=mean(x2),sd_x2=sd(x2))

new <- list()
for (i in 1:nrow(params)) {
  x1 <- rnorm(10,mean=params[i,]$mean_x1,sd=params[i,]$sd_x1)
  x2 <- rnorm(10,mean=params[i,]$mean_x2,sd=params[i,]$sd_x2)
  y <- rep(params[i,]$y,times=10)
  new[[i]] <- data.frame(y,x1,x2)
}

new1 <- bind_rows(new)
write.csv(new1,"rfApp/data/mock_pred.csv",fileEncoding = "utf8",row.names = FALSE)
write.csv(new1[,2:3],"rfApp/data/mock_pred_noY.csv",fileEncoding = "utf8",row.names = FALSE)
ggplot(new1,aes(x=x1,y=x2,color=y)) +
  geom_point()
