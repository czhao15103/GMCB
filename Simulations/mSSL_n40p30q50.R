
#### mSSL, n = 40, p = 30, q = 50 ####

# we use the default choices for hyperparameters and the default ladders for the penalty parameters

library(mSSL)
library(parallel)
library(doParallel)
library(doRNG)

# datasets to use
data <- readRDS("Data_n40p30q50.rds")
y <- data$y.centered
x <- data$x.stand

# simulation settings
p <- ncol(x)
q <- ncol(y[[1]])
n <- nrow(y[[1]])

reps <- length(y) # number of simulation replications

set.seed(325)
start.time <- Sys.time()
ests <- foreach(i = 1:reps, .packages = c("mSSL")) %dorng% {
	out_dpe <- mSSL_dpe(X = x, Y = y[[i]], max_iter = 1e5)
    out_dcpe <- mSSL_dcpe(X = x, Y = y[[i]], max_iter = 1e5)

	savedpe <- paste("mSSL-DPE_rep", i, ".rds", sep = "")
	savedcpe <- paste("mSSL-DCPE_rep", i, ".rds", sep = "")
	saveRDS(out_dpe, file = savedpe)
	saveRDS(out_dcpe, file = savedcpe)
}
end.time <- Sys.time()

time.taken <- end.time - start.time

print(time.taken)




