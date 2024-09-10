
#### GHS n = 40, p = 1, q = 50, precision matrix only ####

library(R.matlab)
library(GMCB)

#### datasets to use ####
data <- readRDS("Data_n40q50_meancovariance_longrange.rds")
true.b <- data$true.b
true.omega <- data$true.omega
true.sigma <- data$true.sigma
reps <- length(data$y)
p <- 1
q <- ncol(true.omega)

#### additional loss functions ####

# scalar quadratic loss
scalarquadloss <- function(b.est, b.actual, omega.actual) {
  b.diff <- b.est - b.actual
  return(sum(diag(b.diff %*% omega.actual %*% t(b.diff))))
}

# stein's loss
kl.loss <- function(omega.est, cov.actual) {
  q <- ncol(omega.est)
  mat.prod <- omega.est %*% cov.actual
  return(sum(diag(mat.prod)) - as.numeric(determinant(mat.prod, logarithm = TRUE)$modulus) - q)
}

#### additional loss functions ####

# stein's loss
kl.loss <- function(omega.est, cov.actual) {
  q <- ncol(omega.est)
  mat.prod <- omega.est %*% cov.actual
  return(sum(diag(mat.prod)) - as.numeric(determinant(mat.prod, logarithm = TRUE)$modulus) - q)
}

#### producing estimates and CIs #### 

# for reading and loading files
filestart <- "n40p1q50_GHS_preconly_rep"

# for saving files
estname <- "GHS_preconly_n40p1q50_ests"

# 2 loss functions each for B and Omega, and 2 estimates each of B and Omega
losses <- matrix(NA, nrow = reps, ncol = 4)
colnames(losses) <- c("F_loss_Omega_postmean",
					  "F_loss_Omega_stein", 
					  "Stein_loss_Omega_postmean", "Stein_loss_Omega_stein")

for (i in 1:reps) {
	cat("i = ", i , "\n")
	filename <- paste(filestart, i, ".mat", sep = "")
	out <- readMat(filename)
	omegas <- lapply(1:dim(out$omega.save)[3], function(x) out$omega.save[,,x])

	### compute estimates
                              
    # Bayes estimator under L_F
    omega.postmean <- Reduce(`+`, omegas)/length(omegas)
                              
    # Bayes estimator under Stein's loss for precision matrix
    cov.mat <- lapply(omegas, solve)
    sigma.postmean <- Reduce(`+`, cov.mat)/length(cov.mat)
    omega.stein <- solve(sigma.postmean)
                              
	estsave <- paste(estname, i, ".rds", sep = "")
	saveRDS(list(omega.postmean = omega.postmean,
                 omega.stein = omega.stein),
			file = estsave)

	### computing the losses

	# squared Frobenius loss, omega post mean
	losses[i,1] <- sum((omega.postmean - true.omega)^2) 

	# squared Frobenius loss, omega stein's 
	losses[i,2] <- sum((omega.stein - true.omega)^2) 

	# Stein's loss, omega post mean
	losses[i,3] <- kl.loss(omega.postmean, true.sigma)

	# Stein's loss, omega stein's 
	losses[i,4] <- kl.loss(omega.stein, true.sigma) 

	### 95% credible intervals
	omega.vech <- vech.list(omegas)
	omega.credible <- sapply(1:ncol(omega.vech), function(x) quantile(omega.vech[,x], probs = c(0.025, 0.975)))

	credsave <- paste("GHS_preconly_n40p1q50_95CredInt", i, ".rds", sep = "")
	saveRDS(omega.credible,
			file = credsave)

	#unlink(filename) # to avoid using too much memory
}

loss_save <- "GHS_preconly_n40p1q50_losses.rds"
saveRDS(losses, file = loss_save)

