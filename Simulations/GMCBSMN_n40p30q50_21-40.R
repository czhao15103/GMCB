
#### GMCB-SMN algorithm simulation, n = 40, p = 30, q = 50 ####

library(GMCB)

# datasets to use
data <- readRDS("Data_n40p30q50.rds")
y <- data$y.centered
x <- data$x.stand

# simulation settings
p <- ncol(x)
q <- ncol(y[[1]])
n <- nrow(y[[1]])

reps <- length(y) # number of simulation replications
iter <- 5e4

# for setting hyperparameters
forscaling <- readRDS("GMCB-SMN_tune14.rds")

# hyperparameters
alpha.prior <- forscaling$priors$alpha
lambda.prior <- forscaling$priors$lambda
tau.prior <- forscaling$priors$tau

rw.scale <- forscaling$rw.scale

clus <- makeCluster(10)
registerDoParallel(clus)

set.seed(943)
ests <- foreach(i = 21:40,
                .packages = c("GMCB"),
                .errorhandling = c("pass")) %dorng% {
	init <- init_freqbridge_nocovariates(y[[i]], alpha.prior, cores = 1)

	priors <- list(lambda = lambda.prior, tau = tau.prior,
			   alpha = alpha.prior, gamma = init$gamma.prior)

	b.init <- solve(crossprod(x), crossprod(x, y[[i]]))
	initial.values <- list(b = as.numeric(matrixcalc::vec(b.init)), 
					   lambda = rep(init$lambda.init[1], p*q),
					   alpha.b = init$alpha.b.init,
					   delta = init$delta.init, 
					   tau = init$tau.init, 
					   gamma = init$gamma.init, 
                       alpha.d = init$alpha.d.init)

     out <- gmcb(y = y[[i]], x = x, 
				 initial.values = initial.values, priors = priors, 
                 rw.scale = rw.scale, iter = iter,
                 algorithm = "SMN")
	 filesave <- paste("run", i, ".rds", sep = "") 
     saveRDS(out, file = filesave)
                            }

stopCluster(clus)






