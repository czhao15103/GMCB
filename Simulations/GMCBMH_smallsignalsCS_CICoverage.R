
#### GMCB-MH algorithm simulation, p = 5, q = 5, n = 100 ####

# This is the dense B, compound symmetric Sigma setting
# X here has an AR(1) matrix 

# the y and x matrices have been standardized already

# functions and packages
library(GMCB)

# datasets to use
data <- readRDS("Data_n100p5q5_smallsignalsCS.rds")
y <- data$y.centered
x <- data$x.stand

# simulation settings
p <- ncol(x)
q <- ncol(y[[1]])
n <- nrow(y[[1]])

reps <- length(y) # number of simulation replications
iter <- 1.5e5 # number of steps for sampler

# hyperparameters
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

clus <- makeCluster(20)
registerDoParallel(clus)

set.seed(363)
ests <- foreach(i = 1:reps,
                .packages = c("GMCB")) %dorng% {
                              out <- gmcb(y = y[[i]], x = x, priors = priors,
                                          rw.scale = rw.scale, iter = iter,
                                          algorithm = "MH")
                              omega <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = TRUE) 
                              
							  b.credible <- sapply(1:ncol(out$mcmc$b), function(x) quantile(out$mcmc$b[,x], probs = c(0.025, 0.975)))

							  omega.vech <- vech.list(omega$prec.mat)
							  omega.credible <- sapply(1:ncol(omega.vech), function(x) quantile(omega.vech[,x], probs = c(0.025, 0.975)))

	credsave <- paste(getwd(), "/", "GMCBMH_smallsignalsCS_CICoverage", i, ".rds", sep = "")
	saveRDS(list(b.credible = b.credible, omega.credible = omega.credible),
			file = credsave)
                            }

stopCluster(clus)



