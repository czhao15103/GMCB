
#### Computation Assessment: n = 100, p = 120, q = 50 ####

library(GMCB)
library(mvtnorm)

n <- 100
p <- 120
q <- 50

nsets <- 50

rho <- 0.7

# true response covariance and precision matrices
y.mats <- matgen(p = q, rho = rho, type = "AR")

# predictor matrix covariance
x.mats <- matgen(p = p, rho = rho, type = "AR")

set.seed(2023)
nonzero <- p*q/2
b <- numeric(p*q)
nonzeroind <- sample(p*q, size = nonzero)
b[nonzeroind] <- runif(nonzero, min = 0.5, max = 2)
b[nonzeroind[1:nonzero/2]] <- -b[nonzeroind[1:nonzero/2]]
true.b <- matrix(b, nrow = p, ncol = q)
write.csv(true.b, file = "b.csv",
          row.names = FALSE)

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
write.csv(x.stand, file = "x.csv",
          row.names = FALSE)

# centered version of y for estimation
y.means <- lapply(y, colMeans)
y.centered <- lapply(1:nsets, function(j) sweep(y[[j]], 2, y.means[[j]]))

values <- list(y = y, x = x, b = true.b, y.centered = y.centered, x.stand = x.stand)
saveRDS(values, file = "Data_n100p120q50_moderatesparsesignalsAR.rds")

# for use with HSGHS
for (i in 1:nsets) {
  # file name 
  y.name <- paste("y", i, ".csv", sep = "")
  
  # save to csv
  write.csv(values$y.centered[[i]], file = y.name, row.names = FALSE)
}
