
#### LOO posterior predictive checks ####

# Note that this code assumes a folder named LOO_GMCBMH_Sitka1988_0112001 exists in the 
# current working directory

data <- readRDS("Sitka1988_arranged.rds")
y.stand <- data$y.stand
x <- data$x

n <- nrow(y.stand)
p <- ncol(x)
q <- ncol(y.stand)

library(GMCB)

#### GMCB-MH ####

###### (0.1, 1, 2, 0.01) ######

forscaling <- readRDS("GMCBMH_Sitka1988_start_0112001.rds")
iter <- 3e5

clus <- makeCluster(20) # replace 20 with number of available cores
registerDoParallel(clus)

set.seed(312)
logpostpred_mh_0112001 <- foreach(i = 1:n, .combine = c, .packages = c("GMCB")) %dorng% {
  # fit the model
  outn <- gmcb(y.stand[-i,], x[-i,], initial.values = forscaling$final.values,
               priors = forscaling$priors, 
              rw.scale = forscaling$rw.scale, iter = iter,
              algorithm = "MH")
  saveRDS(outn, file = paste("LOO_GMCBMH_Sitka1988_0112001/LOO_GMCBMH_Sitka1988_0112001_", i, ".rds", sep = ""))
  
  # posterior predictive distribution values
  bn <- b.to.matrix(outn$mcmc$b, p, q) # list of the B matrices
  tn <- delta.to.matrix(outn$mcmc$delta, q) # list of the T matrices
  tinvn <- lapply(tn, forwardsolve, x = diag(q)) # list of the T^inv matrices
  
  postpredval <- sapply(1:iter, function(z) tinvn[[z]] %*% dnorm(y.stand[i,], mean = t(bn[[z]]) %*% x[i,],
                                                      sd = sqrt(outn$mcmc$gamma[z,])))
  return(log(mean(postpredval)))
}

saveRDS(logpostpred_mh_0112001, file = "LOOPPC_MH_0112001.rds")
