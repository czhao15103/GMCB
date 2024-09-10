
#### n = 40, p = 30, q = 50 Simulations ####

# This file generates additional simulation settings for which n = 40, p = 30, and q = 50

#### packages and settings ####

library(GMCB)
library(mvtnorm)
library(mcmcse)

n <- 40
p <- 30
q <- 50

# number of parameters of interest
p*q + q*(q+1)/2 # 2775

# minimum effective sample size
minESS(p*q + q*(q+1)/2) # 7113  

rho <- 0.7

x.mats <- matgen(p, rho, type = "AR")

#### generating b, omega, x, and y ####

set.seed(14)

# fix x for all datasets
x <- rmvnorm(n, mean = rep(0, p), sigma = x.mats$cov.mat)

# b - 5% non-zero entries, with non-zero entries sampled from U(0.5,2) and sign set randomly
b <- matrix(0, nrow = p, ncol = q)
numnonzero <- p*q * 0.05
magnitudes <- runif(numnonzero, min = 0.5, max = 2)
signs <- sample(c(-1,1), size = numnonzero, replace = TRUE, prob = c(0.5, 0.5))
nonzeroind <- sample(1:(p*q), size = numnonzero)
b[nonzeroind] <- magnitudes * signs

# omega - cliques; 16 groups with 3 members as in Li et al. (2021) setting 2 with q = 50
omega <- diag(q)
ncliques <- 16
cliquesize <- 3
groupmembers <- sample(q, ncliques * cliquesize)
for (i in 1:ncliques) {
  cliquemem <- groupmembers[((i-1)*cliquesize + 1):(cliquesize*i)]
  for (j in 1:cliquesize) {
    for (k in 1:cliquesize) {
      omega[cliquemem[j], cliquemem[k]] <- 0.75
    }
  }
}
diag(omega) <- 1
sigma <- solve(omega)

#### save datasets ####

nsets <- 100

y <- vector("list", length = nsets)

for (i in 1:nsets) {
	y[[i]] <- x %*% b + rmvnorm(n, mean = rep(0, q), sigma = sigma)
}

# standardized version of x to use for estimation
x.means <- colMeans(x)
x.sds <- apply(x, 2, sd)
x.centered <- sweep(x, 2, x.means)
x.stand <- sweep(x.centered, 2, x.sds, FUN = "/")

# centered version of y for estimation
y.means <- lapply(y, colMeans)
y.centered <- lapply(1:nsets, function(j) sweep(y[[j]], 2, y.means[[j]]))

values <- list(y = y, x = x, b = b, y.centered = y.centered, x.stand = x.stand, 
               true.sigma = sigma, true.omega = omega)
saveRDS(values, file = "Data_n40p30q50.rds")

