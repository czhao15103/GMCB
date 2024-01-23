
#### ACF comparison, DP2002 3A XJ2019 1 ####

library(GMCB)
library(mvtnorm)

data <- readRDS("Data_n100p5q5_DP2002_3A_XJ2019_1.rds")

n <- 100
p <- 5
q <- 5

x <- data$x.stand
true.b <- data$b
true.cov <- data$y.cov

# set seed
set.seed(16)

# generate Y
y <- data$x %*% true.b + rmvnorm(n, mean = rep(0, q), sigma = true.cov)

y.means <- colMeans(y)
y <- sweep(y, 2, y.means)

alpha.prior <- c(0.5, 4)
lambda.prior <- c(0.1, 1, 2, 0.01)
tau.prior <- c(1, 1, 40, 0.5)
priors <- list(lambda = lambda.prior, tau = tau.prior, alpha = alpha.prior)

b.scale <- c(rep(0.1, 17), 0.2, rep(0.1, 2), 0.15, 0.3, rep(0.15, 3))
d.scale <- 0.2
alpha.b.scale <- 0.3
alpha.d.scale <- 1
rw.scale <- list(b = b.scale, delta = d.scale, alpha.b = alpha.b.scale,
                 alpha.d = alpha.d.scale)
iter <- 1.5e5 # number of steps for sampler

write.csv(y, file = "acf_data_n100p5q5_DP2002_3A_XJ2019_1_y.csv",
          row.names = FALSE)

outMH <- gmcb(y = y, x = x, priors = priors, rw.scale = rw.scale, iter = iter,
              algorithm = "MH")
saveRDS(outMH, file = "acf_gmcb_mh_n100p5q5_DP2002_3A_XJ2019_1.rds")

# change alpha.b.scale and alpha.d.scale for GMCB-SMN and support of alpha.prior
alpha.prior <- c(0.5, 2)
lambda.prior <- c(0.1, 1, 2, 0.01)
tau.prior <- c(1, 1, 40, 0.5)
priors <- list(lambda = lambda.prior, tau = tau.prior, alpha = alpha.prior)

rw.scale <- list(alpha.b = 0.3, alpha.d = 0.25)

iter <- 2.75e4

outsmn <- gmcb(y = y, x = x, priors = priors, rw.scale = rw.scale, iter = iter,
               algorithm = "SMN")
saveRDS(outsmn, file = "acf_gmcb_smn_n100p5q5_DP2002_3A_XJ2019_1.rds")
