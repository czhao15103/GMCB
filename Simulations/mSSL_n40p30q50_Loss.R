
#### mSSL n = 40, p = 30, q = 50, loss ####

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

#### producing estimates #### 

# for reading and loading files
filestart.dcpe <- "mSSL-DCPE_rep"
filestart.dpe <- "mSSL-DPE_rep"

# 2 loss functions each for B and Omega, and 2 estimates each of B and Omega (DPE and DCPE)
losses <- matrix(NA, nrow = reps, ncol = 8)
colnames(losses) <- c("F_loss_B_DPE", "F_loss_B_DCPE", "F_loss_Omega_DPE",
					  "F_loss_Omega_DCPE", "SQ_loss_B_DPE", "SQ_loss_B_DCPE",
					  "Stein_loss_Omega_DPE", "Stein_loss_Omega_DCPE")

for (i in 1:reps) {
	cat("i = 1", i , "\n")
	filename.dcpe <- paste(filestart.dcpe, i, ".rds", sep = "")
	filename.dpe <- paste(filestart.dpe, i, ".rds", sep = "")
	dcpe <- readRDS(filename.dcpe)
	dpe <- readRDS(filename.dpe)

	### computing the losses

	# squared Frobenius loss, B DPE
	losses[i,1] <- sum((dpe$B - true.b)^2) 

	# squared Frobenius loss, B DCPE
	losses[i,2] <- sum((dcpe$B - true.b)^2) 

	# squared Frobenius loss, omega DPE
	losses[i,3] <- sum((dpe$Omega - true.omega)^2) 

	# squared Frobenius loss, omega DCPE
	losses[i,4] <- sum((dcpe$Omega - true.omega)^2) 

	# scalar quadratic loss, B DPE
	losses[i,5] <- scalarquadloss(dpe$B, true.b, true.omega)

	# scalar quadratic loss, B DCPE
	losses[i,6] <- scalarquadloss(dcpe$B, true.b, true.omega)

	# Stein's loss, omega DPE
	losses[i,7] <- kl.loss(dpe$Omega, true.sigma)

	# Stein's loss, omega DCPE
	losses[i,8] <- kl.loss(dcpe$Omega, true.sigma) 

	unlink(filename.dpe) # to avoid using too much memory
	unlink(filename.dcpe) # to avoid using too much memory
}

loss_save <- "mSSL_n40p30q50_losses.rds"
saveRDS(losses, file = loss_save)

