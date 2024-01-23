
#### Simulation Data, n = 100, p = 5, q = 5 ####

#### Functions, packages, and general common settings ####

library(GMCB)
library(mvtnorm)

p <- 5
q <- 5
n <- 100
x.rho <- 0.7
x.mats <- matgen(p, x.rho, type = "AR")
nsets <- 2000

#### CS, small signals ####

y.mats <- matgen(q, 0.7, type = "CS") 

set.seed(25)

(true.b <- matrix(rnorm(p*q, mean = 2, sd = 0.001), nrow = p))
x <- rmvnorm(n, mean = rep(0, p), sigma = x.mats$cov.mat)

y <- vector("list", length = nsets)

for (i in 1:nsets) {
  y[[i]] <- x %*% true.b + rmvnorm(n, mean = rep(0, q), sigma = y.mats$cov.mat)
}

# standardized version of x to use for estimation
x.means <- colMeans(x)
x.sds <- apply(x, 2, sd)
x.centered <- sweep(x, 2, x.means)
x.stand <- sweep(x.centered, 2, x.sds, FUN = "/")

# centered version of y for estimation
y.means <- lapply(y, colMeans)
y.centered <- lapply(1:nsets, function(j) sweep(y[[j]], 2, y.means[[j]]))

values <- list(y = y, x = x, b = true.b, y.centered = y.centered, x.stand = x.stand)
saveRDS(values, file = "Data_n100p5q5_smallsignalsCS.rds")

#### AR1, moderate, sparse signals ####

y.mats <- matgen(q, 0.7, type = "AR") 

set.seed(26)

true.b <- matrix(rnorm(p*q, mean = 5, sd = 1), nrow = p)
set.to.zero <- sample(1:25, size = 12)
true.b[set.to.zero] <- 0
true.b
x <- rmvnorm(n, mean = rep(0, p), sigma = x.mats$cov.mat)

y <- vector("list", length = nsets)

for (i in 1:nsets) {
  y[[i]] <- x %*% true.b + rmvnorm(n, mean = rep(0, q), sigma = y.mats$cov.mat)
}

# standardized version of x to use for estimation
x.means <- colMeans(x)
x.sds <- apply(x, 2, sd)
x.centered <- sweep(x, 2, x.means)
x.stand <- sweep(x.centered, 2, x.sds, FUN = "/")

# centered version of y for estimation
y.means <- lapply(y, colMeans)
y.centered <- lapply(1:nsets, function(j) sweep(y[[j]], 2, y.means[[j]]))

values <- list(y = y, x = x, b = true.b, y.centered = y.centered, x.stand = x.stand)
saveRDS(values, file = "Data_n100p5q5_moderatesparsesignalsAR.rds")

#### Daniels and Pourahmadi (2002) Scenario 3A, Xiang and Jones (2019) Scenario I ####

# y covariance and precision matrices
gamma <- c(0.5, 0.7, 1, 3, 5)
delta <- matrix(0, nrow = 5, ncol = 5)
for (i in 2:5) {
  for (j in 1:(i-1)) {
    if (j == (i - 1)) {
      delta[i,j] <- -(0.75 + 0.02*j)
    } else if (j == (i - 2)) {
      delta[i,j] <- -0.4
    } else if (j == (i - 3)) {
      delta[i,j] <- -0.2
    } else if (j == (i - 4)) {
      delta[i,j] <- -0.1
    }
  }
}
delta
diag(delta) <- 1
delta
(prec <- t(delta) %*% diag(1/gamma) %*% delta)
(cov.mat <- solve(prec)) # covariance matrix

set.seed(27)

true.b <- matrix(0, nrow = p, ncol = q)
nonzero <- sample(1:25, size = 3)
true.b[nonzero] <- rnorm(3, mean = 15, sd = 3)
true.b
x <- rmvnorm(n, mean = rep(0, p), sigma = x.mats$cov.mat)

y <- vector("list", length = nsets)

for (i in 1:nsets) {
  y[[i]] <- x %*% true.b + rmvnorm(n, mean = rep(0, q), sigma = cov.mat)
}

# standardized version of x to use for estimation
x.means <- colMeans(x)
x.sds <- apply(x, 2, sd)
x.centered <- sweep(x, 2, x.means)
x.stand <- sweep(x.centered, 2, x.sds, FUN = "/")

# centered version of y for estimation
y.means <- lapply(y, colMeans)
y.centered <- lapply(1:nsets, function(j) sweep(y[[j]], 2, y.means[[j]]))

values <- list(y = y, x = x, b = true.b, y.cov = cov.mat, y.prec = prec,
               y.centered = y.centered, x.stand = x.stand)
saveRDS(values, file = "Data_n100p5q5_DP2002_3A_XJ2019_1.rds")
