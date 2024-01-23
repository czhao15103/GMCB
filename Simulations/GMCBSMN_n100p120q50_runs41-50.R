
#### GMCB-SMN Computational Efficiency, n = 100, p = 120, q = 50 ####

# 50 repetitions required more than available time resources - so had to split into two runs
# second set of runs

library(GMCB)

# datasets to use
data <- readRDS("Data_n100p120q50_moderatesparsesignalsAR.rds")
y <- data$y.centered
x <- data$x.stand

# simulation settings
p <- ncol(x)
q <- ncol(y[[1]])
n <- nrow(y[[1]])

reps <- length(y) # number of simulation replications
iter <- 500

# hyperparameters
alpha.prior <- c(0.5, 2)
lambda.prior <- c(0.1, 1, 2, 0.01)
tau.prior <- c(0.1, 1, 2, 0.01) 
priors <- list(lambda = lambda.prior, tau = tau.prior, alpha = alpha.prior)

alpha.b.scale <- 0.15
alpha.d.scale <- 0.25
rw.scale <- list(alpha.b = alpha.b.scale, alpha.d = alpha.d.scale)

#clus <- makeCluster(20)
clus <- makeCluster(10)  # change according to number of available cores
registerDoParallel(clus)

set.seed(100)
ests <- foreach(i = 41:reps,
                .packages = c("GMCB"),
                .errorhandling = c("pass")) %dorng% {
                  out <- gmcb(y = y[[i]], x = x, priors = priors, 
                              rw.scale = rw.scale, iter = iter,
                              algorithm = "SMN")
				  filesave <- paste("GMCBSMN_run", i, 
									".rds", sep = "") 
                  saveRDS(out, file = filesave)
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

stopCluster(clus)

saveRDS(ests, file = "smn_n100p120q50_moderatesparsesignalsAR_41-50.rds")



