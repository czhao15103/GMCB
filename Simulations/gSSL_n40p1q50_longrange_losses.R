
#### gSSL n = 40, p = 1, q = 50 ####

library(GMCB)

#### datasets to use ####
data <- readRDS("Data_n40q50_meancovariance_longrange.rds")
true.omega <- data$true.omega
true.sigma <- data$true.sigma
reps <- length(data$y)
q <- ncol(true.omega)

#### additional loss functions ####

# stein's loss
kl.loss <- function(omega.est, cov.actual) {
  q <- ncol(omega.est)
  mat.prod <- omega.est %*% cov.actual
  return(sum(diag(mat.prod)) - as.numeric(determinant(mat.prod, logarithm = TRUE)$modulus) - q)
}

#### loss #### 

# for reading and loading files
filestart <- "gSSL_preconly_n40p1q50_longrange_rep"

# 2 loss functions for Omega, and 2 estimates of Omega
losses <- matrix(NA, nrow = reps, ncol = 2)
colnames(losses) <- c("F_loss", "Stein_loss")

for (i in 1:reps) {
	cat("i = 1", i , "\n")
	filename <- paste(filestart, i, ".rds", sep = "")
	out <- readRDS(filename)

	### computing the losses

	# squared Frobenius loss, omega post mean
	losses[i,1] <- sum((out$Omega - true.omega)^2) 

	# Stein's loss, omega post mean
	losses[i,2] <- kl.loss(out$Omega, true.sigma)

	unlink(filename) # to avoid using too much memory
}

loss_save <- "gSSL_preconly_n40p1q50_longrange_losses.rds"
saveRDS(losses, file = loss_save)

