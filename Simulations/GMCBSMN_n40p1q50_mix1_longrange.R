
#### n40q50, long range dependence ####

library(GMCB)

data <- readRDS("Data_n40q50_meancovariance_longrange.rds")
y <- data$y
reps <- length(y)
iter <- 5e4 # number of iterations for sampler

# parameter settings
alpha.prior <- c(0.5, 2)
lambda.prior <- c(1, 1, 40, 0.5)
tau.prior <- c(0.1, 1, 2, 0.01)
priors <- list(lambda = lambda.prior, tau = tau.prior,
               alpha = alpha.prior)

rw.scale <- list(alpha.b = 0.2, alpha.d = 0.065)

# for saving files 
filenamestart <- "SMN_n40p1q50_mix1_longrange_rep"

clus <- makeCluster(20)
registerDoParallel(clus)

set.seed(150)
start.time <- Sys.time()
ests <- foreach(i = 1:reps, .packages = c("GMCB")) %dorng% {
	# run algorithm 
	out <- gmcb(y = y[[i]], meanzero = FALSE,  
				priors = priors, rw.scale = rw.scale, iter = iter, algorithm = "SMN")
	filename <- paste(filenamestart, i, ".rds", sep = "")
	saveRDS(out, filename)
}

stopCluster(clus)

