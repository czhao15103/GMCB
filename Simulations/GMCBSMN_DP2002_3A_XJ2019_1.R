
#### GMCB-SMN algorithm simulation, p = 5, q = 5, n = 100 ####

# DP 2002 3A XJ 2019 1

library(GMCB)

# datasets to use
data <- readRDS("Data_n100p5q5_DP2002_3A_XJ2019_1.rds")
y <- data$y.centered
x <- data$x.stand

# simulation settings
p <- ncol(x)
q <- ncol(y[[1]])
n <- nrow(y[[1]])

reps <- length(y) # number of simulation replications
iter <- 2.75e4 # number of steps for sampler

# hyperparameters
alpha.prior <- c(0.5, 2)
lambda.prior <- c(0.1, 1, 2, 0.01)
tau.prior <- c(1, 1, 40, 0.5)
priors <- list(lambda = lambda.prior, tau = tau.prior, alpha = alpha.prior)

rw.scale <- list(alpha.b = 0.3, alpha.d = 0.25)

clus <- makeCluster(20)
registerDoParallel(clus)

set.seed(0)
start.time <- Sys.time()
ests <- foreach(i = 1:reps,
                .packages = c("GMCB"),
                .errorhandling = c("pass")) %dorng% {
                  out <- gmcb(y = y[[i]], x = x, priors = priors, 
                              rw.scale = rw.scale, iter = iter,
                              algorithm = "SMN")
                  omega <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = TRUE) 
                  
                  # Bayes estimator under L_F
                  omega.postmean <- Reduce(`+`, omega$prec.mat)/length(omega$prec.mat)
                  b.postmean.vec <- colMeans(out$mcmc$b)
                  b.postmean <- matrix(b.postmean.vec, nrow = p, ncol = q)
                  
                  # Bayes estimator under Stein's loss for precision matrix
                  sigma.postmean <- Reduce(`+`, omega$cov.mat)/length(omega$cov.mat)
                  omega.stein <- solve(sigma.postmean, tol = 1e-30)
                  
                  # Bayes estimator under scalar quadratic loss for B
                  b.matrix <- b.to.matrix(b = out$mcmc$b, p = p, q = q)
                  b.times.sigmainv <- lapply(1:length(omega$prec.mat), 
                                             function(x) b.matrix[[x]] %*% omega$prec.mat[[x]])
                  b.times.sigmainv.postmean <- Reduce(`+`, b.times.sigmainv)/length(b.times.sigmainv)
                  b.quad <- b.times.sigmainv.postmean %*% solve(omega.postmean, tol = 1e-30)
                              
                              return(list(b.postmean = b.postmean, omega.postmean = omega.postmean,
                                          b.quad = b.quad, omega.stein = omega.stein))
                            }
end.time <- Sys.time()

stopCluster(clus)

time.taken <- end.time - start.time

out <- list(ests = ests, time = time.taken)

saveRDS(out, file = "gmcb_smn_n100p5q5_DP2002_3A_XJ2019_1.rds")




