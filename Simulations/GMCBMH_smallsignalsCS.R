
#### GMCB-MH algorithm simulation, p = 5, q = 5, n = 100 ####

# This is the dense B, compound symmetric Sigma setting

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
start.time <- Sys.time()
ests <- foreach(i = 1:reps,
                .packages = c("GMCB")) %dorng% {
                              out <- gmcb(y = y[[i]], x = x, priors = priors,
                                          rw.scale = rw.scale, iter = iter,
                                          algorithm = "MH")
                              omega <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = TRUE) 
                              
                              # Bayes estimator under L_F
                              omega.postmean <- Reduce(`+`, omega$prec.mat)/length(omega$prec.mat)
                              b.postmean.vec <- colMeans(out$mcmc$b)
                              b.postmean <- matrix(b.postmean.vec, nrow = p, ncol = q)
                              
                              # Bayes estimator under Stein's loss for precision matrix
                              sigma.postmean <- Reduce(`+`, omega$cov.mat)/length(omega$cov.mat)
                              omega.stein <- solve(sigma.postmean)
                              
                              # Bayes estimator under scalar quadratic loss for B
                              b.matrix <- b.to.matrix(b = out$mcmc$b, p = p, q = q)
                              b.times.sigmainv <- lapply(1:length(omega$prec.mat), 
                                                         function(x) b.matrix[[x]] %*% omega$prec.mat[[x]])
                              b.times.sigmainv.postmean <- Reduce(`+`, b.times.sigmainv)/length(b.times.sigmainv)
                              b.quad <- b.times.sigmainv.postmean %*% solve(omega.postmean)
                              
                              return(list(b.postmean = b.postmean, omega.postmean = omega.postmean,
                                          b.quad = b.quad, omega.stein = omega.stein))
                            }
end.time <- Sys.time()

stopCluster(clus)

time.taken <- end.time - start.time

out <- list(ests = ests, time = time.taken)

saveRDS(out, file = "gmcb_mh_n100p5q5_smallsignalsCS.rds")



