
#### Mean-covariance estimation ####

library(GMCB)
library(mvtnorm)
library(mcmcse)

n <- 40
q <- 50

#### long range dependence ####

h <- 0.7
sigma <- matrix(nrow = q, ncol = q)
for (i in 1:q) {
  for (j in 1:q) {
    sigma[i,j] <- 0.5 * (abs((abs(i - j) + 1))^(2*h) - 2*abs(i - j)^(2*h) + abs(abs(i - j) - 1)^(2*h))
  }
}

omega <- solve(sigma)

set.seed(28)

# mean 
mean.mats <- matgen(q, 0.5, type = "CS")
mu <- rmvnorm(1, mean = rep(0, q), sigma = mean.mats$cov.mat)
mu <- sort(as.numeric(mu))

nsets <- 100

y <- vector("list", length = nsets)

for (i in 1:nsets) {
  y[[i]] <- rmvnorm(n, mean = mu, sigma = sigma)
}

values <- list(y = y, true.b = mu, true.sigma = sigma, true.omega = omega)
saveRDS(values, 
        file = "Data_n40q50_meancovariance_longrange.rds")

