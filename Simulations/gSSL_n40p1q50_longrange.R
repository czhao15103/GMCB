
#### gSSL prec only, n = 40, p = 1, q = 50 ####

# we use the default choices for hyperparameters and the default ladders for the penalty parameters

library(mSSL)
library(parallel)
library(doParallel)
library(doRNG)

# datasets to use
data <- readRDS("Data_n40q50_meancovariance_longrange.rds")
y.means <- lapply(data$y, colMeans)
y <- lapply(1:length(data$y), function(j) sweep(data$y[[j]], 2, y.means[[j]]))

# simulation settings
q <- ncol(y[[1]])
n <- nrow(y[[1]])

reps <- length(y) # number of simulation replications

clus <- makeCluster(20)
registerDoParallel(clus)

set.seed(33)
ests <- foreach(i = 1:reps, .packages = c("mSSL")) %dorng% {
	out <- gSSL(Y = y[[i]], eps = 1e-3, max_iter = 1e5)

	savefile <- paste("gSSL_preconly_n40p1q50_longrange_rep", i, ".rds", sep = "")
	saveRDS(out, file = savefile)
}

stopCluster(clus)

