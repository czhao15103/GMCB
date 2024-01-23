
#### DP2002_3A_XJ2019_1 ####

# we use the default choices for hyperparameters and the default ladders for the penalty parameters

library(mSSL)

library(parallel)
library(doParallel)
library(doRNG)

# datasets to use
data <- readRDS("Data_n100p5q5_DP2002_3A_XJ2019_1.rds")
y <- data$y.centered
x <- data$x.stand

# simulation settings
p <- ncol(x)
q <- ncol(y[[1]])
n <- nrow(y[[1]])

reps <- length(y) # number of simulation replications

clus <- makeCluster(20)
registerDoParallel(clus)

set.seed(4)
start.time <- Sys.time()
ests <- foreach(i = 1:reps,
                .packages = c("mSSL")) %dorng% {
                  out_dpe <- mSSL_dpe(X = x, Y = y[[i]])
                  out_dcpe <- mSSL_dcpe(X = x, Y = y[[i]])
                  
                  return(list(b.dpe = out_dpe$B, b.dcpe = out_dcpe$B,
                              omega.dpe = out_dpe$Omega, omega.dcpe = out_dcpe$Omega))
                            }
end.time <- Sys.time()

stopCluster(clus)

time.taken <- end.time - start.time

out <- list(ests = ests, time = time.taken)

saveRDS(out, file = "mSSL_n100p5q5_DP2002_3A_XJ2019_1.rds")
