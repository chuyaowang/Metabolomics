cat("\014")
rm(list = ls())

library(MASS)
library(mixOmics)

# Generate fake dataset using 2 multivariate normal distributions
set.seed(666) # For consistency
numvar <- 2000 # Simulate data 2000 variables
samplesize <- 100 # Total sample size

mean1 <- mvrnorm(n = 1/2*samplesize,
                 mu = )
