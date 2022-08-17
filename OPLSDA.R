cat("\014")
rm(list = ls())

library(ropls)
library(dplyr)
filename <- "D:\\Files\\SH\\Metabo\\data\\sample.csv"
sampledata <- read.csv(filename,
                       header=TRUE)

# Optional removing of the first outlier sample
sampledata <- sampledata[2:6,]

# Get names and group name for convenience
samplename <- sampledata$Sample.Name
samplegroup <- sampledata$Group

# Get only the numerical part
sampleprocessed <- sampledata[,3:length(sampledata)]

# Optional horizontal normalization
sampleprocessed <- apply(sampleprocessed,1,function(x) {
  . <- norm(x,type="2")
  return(x/.)
}) %>%
  t

oplsresult <- opls(x = sampleprocessed,
                   y = samplegroup %>% as.factor,
                   predI = 1,
                   orthoI = NA,
                   crossvalI = 5)
