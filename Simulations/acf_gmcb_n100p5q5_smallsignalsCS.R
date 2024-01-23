
#### ACF comparison, small signals, CS ####

library(GMCB)
library(mvtnorm)

data <- readRDS("Data_n100p5q5_smallsignalsCS.rds")

n <- 100
p <- 5
q <- 5

x <- data$x.stand
true.b <- data$b

# set seed
set.seed(32)

rho <- 0.7
y.mats <- matgen(q, rho, sigma = 1, type = "CS")

# generate Y
y <- data$x %*% true.b + rmvnorm(n, mean = rep(0, q), sigma = y.mats$cov.mat)

y.means <- colMeans(y)
y <- sweep(y, 2, y.means)

alpha.prior <- c(0.5, 4)
lambda.prior <- c(1, 1, 40, 0.5)
tau.prior <- c(1, 1, 40, 0.5)
priors <- list(lambda = lambda.prior, tau = tau.prior, alpha = alpha.prior)

b.scale <- 0.15
d.scale <- 0.1
alpha.b.scale <- 0.3
alpha.d.scale <- 0.7
rw.scale <- list(b = b.scale, delta = d.scale, alpha.b = alpha.b.scale,
                 alpha.d = alpha.d.scale)
iter <- 1.5e5 # number of steps for sampler

write.csv(y, file = "acf_data_n100p5q5_smallsignalsCS_y.csv",
          row.names = FALSE)

outMH <- gmcb(y = y, x = x, priors = priors, rw.scale = rw.scale, iter = iter,
              algorithm = "MH")
saveRDS(outMH, file = "acf_gmcb_mh_n100p5q5_smallsignalsCS.rds")

# change alpha.b.scale and alpha.d.scale for GMCB-SMN and support of alpha.prior
alpha.prior <- c(0.5, 2)
lambda.prior <- c(1, 1, 40, 0.5)
tau.prior <- c(1, 1, 40, 0.5)
priors <- list(lambda = lambda.prior, tau = tau.prior, alpha = alpha.prior)

alpha.b.scale <- 0.3
alpha.d.scale <- 0.3
rw.scale <- list(alpha.b = alpha.b.scale, alpha.d = alpha.d.scale)

iter <- 2.5e4

outsmn <- gmcb(y = y, x = x, priors = priors, rw.scale = rw.scale, iter = iter,
               algorithm = "SMN")
saveRDS(outsmn, file = "acf_gmcb_smn_n100p5q5_smallsignalsCS.rds")
