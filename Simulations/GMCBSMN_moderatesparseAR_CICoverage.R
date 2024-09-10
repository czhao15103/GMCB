
#### GMCB-SMN algorithm simulation, p = 5, q = 5, n = 100 ####

# X here has an AR(1) matrix and we standardize it

library(GMCB)

# datasets to use
data <- readRDS("Data_n100p5q5_moderatesparsesignalsAR.rds")
y <- data$y.centered
x <- data$x.stand

# simulation settings
p <- ncol(x)
q <- ncol(y[[1]])
n <- nrow(y[[1]])

reps <- length(y) # number of simulation replications
iter <- 2.5e4 # number of steps for sampler

# hyperparameters
alpha.prior <- c(0.5, 2)
lambda.prior <- c(0.1, 1, 2, 0.01)
tau.prior <- c(0.1, 1, 2, 0.01) 
priors <- list(lambda = lambda.prior, tau = tau.prior, alpha = alpha.prior)

alpha.b.scale <- 0.15
alpha.d.scale <- 0.25
rw.scale <- list(alpha.b = alpha.b.scale, alpha.d = alpha.d.scale)

clus <- makeCluster(20)
registerDoParallel(clus)

set.seed(91)
ests <- foreach(i = 1:reps,
                .packages = c("GMCB"),
                .errorhandling = c("pass")) %dorng% {
                  out <- gmcb(y = y[[i]], x = x, priors = priors, 
                              rw.scale = rw.scale, iter = iter,
                              algorithm = "SMN")
                  omega <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = TRUE) 
                  
                  							  b.credible <- sapply(1:ncol(out$mcmc$b), function(x) quantile(out$mcmc$b[,x], probs = c(0.025, 0.975)))

							  omega.vech <- vech.list(omega$prec.mat)
							  omega.credible <- sapply(1:ncol(omega.vech), function(x) quantile(omega.vech[,x], probs = c(0.025, 0.975)))

	credsave <- paste(getwd(), "/", "GMCBSMN_moderatesparseAR_CICoverage", i, ".rds", sep = "")
	saveRDS(list(b.credible = b.credible, omega.credible = omega.credible),
			file = credsave)
                            }

stopCluster(clus)



