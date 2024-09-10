
#### HSGHS n = 40, p = 30, q = 50 ####

library(R.matlab)
library(GMCB)

#### datasets to use ####
data <- readRDS("Data_n40p30q50.rds")
true.b <- data$b
true.omega <- data$true.omega
true.sigma <- data$true.sigma
reps <- length(data$y)
p <- nrow(true.b)
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

#### producing estimates and CIs #### 

# for reading and loading files
filestart <- "n40p30q50_HSGHS"

# for saving files
estname <- "HSGHS_n40p30q50_ests"

# 2 loss functions each for B and Omega, and 2 estimates each of B and Omega
losses <- matrix(NA, nrow = reps, ncol = 8)
colnames(losses) <- c("F_loss_B_postmean", "F_loss_B_quad", "F_loss_Omega_postmean",
					  "F_loss_Omega_stein", "SQ_loss_B_postmean", "SQ_loss_B_quad",
					  "Stein_loss_Omega_postmean", "Stein_loss_Omega_stein")

for (i in 1:reps) {
	cat("i = ", i , "\n")
	filename <- paste(filestart, i, ".mat", sep = "")
	out <- readMat(filename)
	omegas <- lapply(1:dim(out$omega.save)[3], function(x) out$omega.save[,,x])

	### compute estimates
                              
    # Bayes estimator under L_F
    omega.postmean <- Reduce(`+`, omegas)/length(omegas)
    b.postmean <- apply(out$beta.save, c(1,2), mean)
                              
    # Bayes estimator under Stein's loss for precision matrix
    cov.mat <- lapply(omegas, solve)
    sigma.postmean <- Reduce(`+`, cov.mat)/length(cov.mat)
    omega.stein <- solve(sigma.postmean)
                              
    # Bayes estimator under scalar quadratic loss for B
    b.times.sigmainv <- lapply(1:length(omegas), 
                               function(x) out$beta.save[,,x] %*% omegas[[x]])
    b.times.sigmainv.postmean <- Reduce(`+`, b.times.sigmainv)/length(b.times.sigmainv)
    b.quad <- b.times.sigmainv.postmean %*% solve(omega.postmean)

	estsave <- paste(estname, i, ".rds", sep = "")
	saveRDS(list(b.postmean = b.postmean, omega.postmean = omega.postmean,
                 b.quad = b.quad, omega.stein = omega.stein),
			file = estsave)

	### computing the losses

	# squared Frobenius loss, B post mean
	losses[i,1] <- sum((b.postmean - true.b)^2) 

	# squared Frobenius loss, B quad loss
	losses[i,2] <- sum((b.quad - true.b)^2) 

	# squared Frobenius loss, omega post mean
	losses[i,3] <- sum((omega.postmean - true.omega)^2) 

	# squared Frobenius loss, omega stein's 
	losses[i,4] <- sum((omega.stein - true.omega)^2) 

	# scalar quadratic loss, B post mean
	losses[i,5] <- scalarquadloss(b.postmean, true.b, true.omega)

	# scalar quadratic loss, B quad loss
	losses[i,6] <- scalarquadloss(b.quad, true.b, true.omega)

	# Stein's loss, omega post mean
	losses[i,7] <- kl.loss(omega.postmean, true.sigma)

	# Stein's loss, omega stein's 
	losses[i,8] <- kl.loss(omega.stein, true.sigma) 

	### 95% credible intervals
	b.vec <- t(sapply(1:dim(out$beta.save)[3], function(x) matrixcalc::vec(out$beta.save[,,x])))
	b.credible <- sapply(1:ncol(b.vec), function(x) quantile(b.vec[,x], probs = c(0.025, 0.975)))

	omega.vech <- vech.list(omegas)
	omega.credible <- sapply(1:ncol(omega.vech), function(x) quantile(omega.vech[,x], probs = c(0.025, 0.975)))

	credsave <- paste("HSGHS95CredInt", i, ".rds", sep = "")
	saveRDS(list(b.credible = b.credible, omega.credible = omega.credible),
			file = credsave)

	unlink(filename) # to avoid using too much memory
}

loss_save <- "HSGHS_n40p30q50_losses.rds"
saveRDS(losses, file = loss_save)


