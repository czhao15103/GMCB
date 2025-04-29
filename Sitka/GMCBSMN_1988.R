
#### GMCB-SMN: Sitka, 1988 ####

data <- readRDS("Sitka1988_arranged.rds")
y.stand <- data$y.stand
x <- data$x

n <- nrow(y.stand)
p <- ncol(x)
q <- ncol(y.stand)

library(mcmcse) 
minESS(p*q + q*(q+1)/2) # 8637 

library(GMCB)

#### (1,1,40,0.5) ####

out <- readRDS("GMCBSMN_Sitka1988_start_114005.rds")

out <- gmcb(out, iter = 1e5)

out$acceptances$alpha_b_accept_rate # 0.36738
par(mfrow = c(1,2))
plot(ts(out$mcmc$alpha_b[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b)]), main = "alpha.b") 
acf(out$mcmc$alpha_b, lag.max = 100) 

out$acceptances$alpha_d_accept_rate # 0.36236
par(mfrow = c(1,2))
plot(ts(out$mcmc$alpha_d[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b)]), main = "alpha.d")
acf(out$mcmc$alpha_d, lag.max = 100) 

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$b)) {
  plot(ts(out$mcmc$b[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b), i]), main = paste("b", i, sep = ""))
  acf(out$mcmc$b[,i], lag.max = 100)
}

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$lambda)) {
  plot(ts(out$mcmc$lambda[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b),i]), main = paste("lambda", i, sep = ""))
  acf(out$mcmc$lambda[-1,i], lag.max = 100)
}

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$delta)) {
  plot(ts(out$mcmc$delta[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b), i]), main = paste("delta", i, sep = ""))
  acf(out$mcmc$delta[,i], lag.max = 100)
}

par(mfrow = c(1,2)) 
for (i in 1:ncol(out$mcmc$tau)) {
  plot(ts(out$mcmc$tau[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b),i]), main = paste("tau", i, sep = ""))
  acf(out$mcmc$tau[-1,i], lag.max = 100)
} 

par(mfrow = c(1,2)) 
for (i in 1:q) {
  plot(ts(out$mcmc$gamma[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b),i]), main = paste("gamma", i, sep = " "))
  acf(out$mcmc$gamma[,i], lag.max = 100)
}

omega <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = FALSE)
omega.mat <- vech.list(omega)
output.of.interest <- cbind(out$mcmc$b, omega.mat)
multiESS(output.of.interest, r = 1) # 48426.22

cltcovest <- mcse.multi(output.of.interest, r = 1)
cltcovest.cov <- cltcovest$cov
sum(diag(cltcovest.cov)) # 693.7166
eigen(cltcovest.cov, only.values = TRUE)$values 
# largest = 4.509285e+02, smallest = 7.271002e-05

steinquad <- ess_steinquad(out$mcmc$b, out$mcmc$delta, out$mcmc$gamma) 
steinquad$multiESS # 43078.55
sum(diag(steinquad$asymp.cov)) # 606.3975
eigen(steinquad$asymp.cov, only.values = TRUE)$values 
# largest =  3.689288e+02, smallest = 6.483410e-05

saveRDS(out, file = "GMCBSMN_Sitka1988_output_114005.rds")

#### (0.1, 1, 2, 0.01) ####

out <- readRDS("GMCBSMN_Sitka1988_start_0112001.rds")

out <- gmcb(out, iter = 1e5)

out$acceptances$alpha_b_accept_rate # 0.41648
par(mfrow = c(1,2))
plot(ts(out$mcmc$alpha_b[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b)]), main = "alpha.b") 
acf(out$mcmc$alpha_b, lag.max = 100) 

out$acceptances$alpha_d_accept_rate # 0.52478
par(mfrow = c(1,2))
plot(ts(out$mcmc$alpha_d[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b)]), main = "alpha.d")
acf(out$mcmc$alpha_d, lag.max = 100) 

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$b)) {
  plot(ts(out$mcmc$b[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b), i]), main = paste("b", i, sep = ""))
  acf(out$mcmc$b[,i], lag.max = 100)
}

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$lambda)) {
  plot(ts(out$mcmc$lambda[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b),i]), main = paste("lambda", i, sep = ""))
  acf(out$mcmc$lambda[-1,i], lag.max = 100)
}

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$delta)) {
  plot(ts(out$mcmc$delta[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b), i]), main = paste("delta", i, sep = ""))
  acf(out$mcmc$delta[,i], lag.max = 100)
}

par(mfrow = c(1,2)) 
for (i in 1:ncol(out$mcmc$tau)) {
  plot(ts(out$mcmc$tau[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b),i]), main = paste("tau", i, sep = ""))
  acf(out$mcmc$tau[-1,i], lag.max = 100)
} 

par(mfrow = c(1,2)) 
for (i in 1:q) {
  plot(ts(out$mcmc$gamma[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b),i]), main = paste("gamma", i, sep = " "))
  acf(out$mcmc$gamma[,i], lag.max = 100)
}

omega <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = FALSE)
omega.mat <- vech.list(omega)
output.of.interest <- cbind(out$mcmc$b, omega.mat)
multiESS(output.of.interest, r = 1) # 49945.69

cltcovest <- mcse.multi(output.of.interest, r = 1)
cltcovest.cov <- cltcovest$cov
sum(diag(cltcovest.cov)) # 732.5979
eigen(cltcovest.cov, only.values = TRUE)$values 
# largest = 5.476973e+02, smallest = 2.234489e-05

steinquad <- ess_steinquad(out$mcmc$b, out$mcmc$delta, out$mcmc$gamma) 
steinquad$multiESS # 41517.7
sum(diag(steinquad$asymp.cov)) # 586.0613
eigen(steinquad$asymp.cov, only.values = TRUE)$values 
# largest =  4.095518e+02, smallest = 2.403155e-05

saveRDS(out, file = "GMCBSMN_Sitka1988_output_0112001.rds")
