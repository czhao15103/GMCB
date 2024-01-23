
#### GMCB-MH Computational Efficiency, p = 5, q = 5, n = 100 ####

# DP 2002 3A XJ 2019 1

# the y and x matrices have been standardized already

# functions and packages
library(GMCB)
library(mcmcse)

# datasets to use
data <- readRDS("Data_n100p5q5_DP2002_3A_XJ2019_1.rds")
y <- data$y.centered
x <- data$x.stand

# simulation settings
p <- ncol(x)
q <- ncol(y[[1]])
n <- nrow(y[[1]])

reps <- length(y) # number of simulation replications
iter <- 1e5 # number of steps for sampler

# hyperparameters
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

clus <- makeCluster(20)
registerDoParallel(clus)

set.seed(111)
start.time <- Sys.time()
ests <- foreach(i = 1:reps,
                .packages = c("GMCB", "mcmcse")) %dorng% {
                  out <- gmcb(y = y[[i]], x = x, priors = priors,
                              rw.scale = rw.scale, iter = iter,
                              algorithm = "MH")
                  
                  # computing multivariate ESS for Bayes estimate under L_F
                  omega <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = FALSE) 
                  
				  omega.mat <- vech.list(omega)
				  output.of.interest <- cbind(out$b, omega.mat)
				  ess_postmean <- multiESS(output.of.interest, r = 1) 

				  cltcovest.cov <- mcse.multi(output.of.interest, r = 1)$cov
				  cltcovest_eigens <- eigen(cltcovest.cov, only.values = TRUE)$values 
				  ess_postmean_cltcovest_largesteigen <- cltcovest_eigens[1]
				  ess_postmean_cltcovest_smallesteigen <- cltcovest_eigens[length(cltcovest_eigens)]
				  
				  # computing multivariate ESS for Bayes estimate under L_Q + L_S
				  steinquadess <- ess_steinquad(out$mcmc$b, out$mcmc$delta, out$mcmc$gamma)
				  steinquad_ceigens <- eigen(steinquadess$asymp.cov, only.values = TRUE)$values 
				  ess_steinquad_cltcovest_largesteigen <- steinquad_ceigens[1]
				  ess_steinquad_cltcovest_smallesteigen <- steinquad_ceigens[length(steinquad_ceigens)]
                              
                  return(list(postmean_ess = ess_postmean, postmean_covtr = sum(cltcovest_eigens),
                              postmean_largesteigen = ess_postmean_cltcovest_largesteigen,
                              postmean_smallesteigen = ess_postmean_cltcovest_smallesteigen,
                              steinquad_ess = steinquadess$multiESS, 
                              steinquad_covtr = sum(steinquad_ceigens) ,
                              steinquad_largesteigen = ess_steinquad_cltcovest_largesteigen,
                              steinquad_smallesteigen = ess_steinquad_cltcovest_smallesteigen,
										          time = out$timing))
  							  }
end.time <- Sys.time()

stopCluster(clus)

time.taken <- end.time - start.time

out <- list(ests = ests, time = time.taken)

saveRDS(out, file = "gmcb_mh_n100p5q5_DP2002_3A_XJ2019_1_ComputationAssessment.rds")



