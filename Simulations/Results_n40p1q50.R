
#### n = 40, p = 1, q = 50 ####

#### sample mean ####

data <- readRDS("Data_n40q50_meancovariance_longrange.rds")
reps <- length(data$y)
p <- 1
q <- ncol(data$y[[1]])

scalarquadloss <- function(b.est, b.actual, omega.actual) {
  b.diff <- b.est - b.actual
  return(sum(diag(b.diff %*% omega.actual %*% t(b.diff))))
}

sampmean <- matrix(nrow = reps, ncol = q)
for (i in 1:reps) {
	sampmean[i,] <- colMeans(data$y[[i]])
}

sampmeanloss <- matrix(nrow = reps, ncol = 2)
colnames(sampmeanloss) <- c("F_loss", "SQ_loss")
for (i in 1:reps) {
	sampmeanloss[i,1] <- sum((sampmean[i,] - data$true.b)^2)
	sampmeanloss[i,2] <- scalarquadloss(matrix(sampmean[i,], nrow = p, ncol = q),
										data$true.b, data$true.omega)
}


# losses
smn_mix1 <- readRDS("SMN_n40p1q50_mix1_longrange_losses.rds")
hsghs <- readRDS("HSGHS_n40p1q50_losses.rds")
ghs <- readRDS("GHS_preconly_n40p1q50_losses.rds")
gssl <- readRDS("gSSL_preconly_n40p1q50_longrange_losses.rds")

rbind(smn_mix1 = colMeans(smn_mix1), 	  hsghs = colMeans(hsghs))

colMeans(ghs)

colMeans(gssl)

colMeans(sampmeanloss)


allsds <- c(sapply(1:4, function(x) sd(smn_mix1[,x])/nrow(smn_mix1)),
	  sapply(1:4, function(x) sd(hsghs[,x])/nrow(hsghs)),
	  sapply(1:2, function(x) sd(ghs[,x])/nrow(ghs)), sd(gssl[,1])/nrow(gssl))
max(allsds)
