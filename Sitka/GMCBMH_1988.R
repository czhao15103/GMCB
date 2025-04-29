
#### GMCB-MH: Sitka, 1988 ####

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

out <- readRDS("GMCBMH_Sitka1988_start_114005.rds")

out <- gmcb(out, iter = 3e5)
saveRDS(out, file = "GMCBMH_Sitka1988_output_114005.rds")

out$acceptances$alpha_b_accept_rate # 0.3325833
par(mfrow = c(1,2))
plot(ts(out$mcmc$alpha_b[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b)]), main = "alpha.b") 
acf(out$mcmc$alpha_b, lag.max = 100) 

out$acceptances$alpha_d_accept_rate # 0.4868767
par(mfrow = c(1,2))
plot(ts(out$mcmc$alpha_d[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b)]), main = "alpha.d")
acf(out$mcmc$alpha_d, lag.max = 100) 

out$acceptances$b_accept_rate
# 0.3240233 0.3541333 0.3799033 0.4381133 0.4099600 0.4331900 0.3869733
# 0.4124300 0.4504067 0.5435367

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$b)) {
  plot(ts(out$mcmc$b[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b), i]), main = paste("b", i, sep = ""))
  acf(out$mcmc$b[,i], lag.max = 100, main = paste("accept = ", round(out$acceptances$b_accept_rate[i], digits = 2)))
}

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$lambda)) {
  plot(ts(out$mcmc$lambda[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b),i]), main = paste("lambda", i, sep = ""))
  acf(out$mcmc$lambda[-1,i], lag.max = 100)
}

out$acceptances$delta_accept_rate
# 0.4701133 0.3405833 0.3230633 0.5711067 0.4022400 0.5054433 0.6050433 0.4894200 0.5012900 0.3823133

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$delta)) {
  plot(ts(out$mcmc$delta[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b), i]), main = paste("delta", i, sep = ""))
  acf(out$mcmc$delta[,i], lag.max = 100, main = paste("accept = ", round(out$acceptances$delta_accept_rate[i], digits = 2)))
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
multiESS(output.of.interest, r = 1) # 10871.9

cltcovest <- mcse.multi(output.of.interest, r = 1)
cltcovest.cov <- cltcovest$cov
sum(diag(cltcovest.cov)) # 10241.66
eigen(cltcovest.cov, only.values = TRUE)$values 
# largest = 5.219296e+03, smallest = 1.274851e-03

steinquad <- ess_steinquad(out$mcmc$b, out$mcmc$delta, out$mcmc$gamma) 
steinquad$multiESS # 10425.95
sum(diag(steinquad$asymp.cov)) # 9708.787
eigen(steinquad$asymp.cov, only.values = TRUE)$values 
# largest =  4.270222e+03, smallest = 1.113847e-03

#### (0.1, 1, 2, 0.01) ####

out <- readRDS("GMCBMH_Sitka1988_start_0112001.rds")

out <- gmcb(out, iter = 3e5)

out$acceptances$b_accept_rate 
# 0.3863933 0.4171000 0.3275600 0.3307633 0.2941467 0.2998967 0.3132633
# 0.3254733 0.2741033 0.3339567

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$b)) {
  plot(ts(out$mcmc$b[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b), i]), main = paste("b", i, sep = ""))
  acf(out$mcmc$b[,i], lag.max = 100,
      main = paste("accept = ", out$acceptances$b_accept_rate[i]))
}

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$lambda)) {
  plot(ts(out$mcmc$lambda[(nrow(out$mcmc$lambda) - 999):nrow(out$mcmc$lambda), i]), main = paste("lambda", i, sep = ""))
  acf(out$mcmc$lambda[-1,i], lag.max = 100)
}

out$acceptances$delta_accept_rate
# 0.4716367 0.4821633 0.5289867 0.4033500 0.4032567 0.5621700 0.3925300
# 0.3914833 0.3928867 0.4366500

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$delta)) {
  plot(ts(out$mcmc$delta[(nrow(out$mcmc$delta) - 999):nrow(out$mcmc$delta), i]), main = paste("delta", i, sep = ""))
  acf(out$mcmc$delta[,i], lag.max = 100, main = paste("accept = ", out$acceptances$delta_accept_rate[i]))
}

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$tau)) {
  plot(ts(out$mcmc$tau[(nrow(out$mcmc$tau) - 999):nrow(out$mcmc$tau), i]), main = paste("tau", i, sep = ""))
  acf(out$mcmc$tau[-1,i], lag.max = 100)
}

par(mfrow = c(1,2))
for (i in 1:ncol(out$mcmc$gamma)) {
  plot(ts(out$mcmc$gamma[(nrow(out$mcmc$gamma) - 999):nrow(out$mcmc$gamma), i]), main = paste("gamma", i, sep = ""))
  acf(out$mcmc$gamma[,i], lag.max = 100)
}

out$acceptances$alpha_b_accept_rate # 0.4035567
par(mfrow = c(1,2))
plot(ts(out$mcmc$alpha_b[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b)]), main = "alpha.b") 
acf(out$mcmc$alpha_b, lag.max = 100, 
    main = paste("accept = ", out$acceptances$alpha_b_accept_rate)) # better

out$acceptances$alpha_d_accept_rate # 0.6857967
par(mfrow = c(1,2))
plot(ts(out$mcmc$alpha_d[(nrow(out$mcmc$b) - 999):nrow(out$mcmc$b)]), main = "alpha.d")
acf(out$mcmc$alpha_d, lag.max = 100, 
    main = paste("accept = ", out$acceptances$alpha_d_accept_rate))# pretty good

omega <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = FALSE)
omega.mat <- vech.list(omega)
output.of.interest <- cbind(out$mcmc$b, omega.mat)
multiESS(output.of.interest, r = 1) # 13749.51

cltcovest <- mcse.multi(output.of.interest, r = 1)
cltcovest.cov <- cltcovest$cov
sum(diag(cltcovest.cov)) # 7825.772
eigen(cltcovest.cov, only.values = TRUE)$values 
# largest = 4.396771e+03, smallest = 5.107978e-04

steinquad <- ess_steinquad(out$mcmc$b, out$mcmc$delta, out$mcmc$gamma) 
steinquad$multiESS # 13450.34
sum(diag(steinquad$asymp.cov)) # 6574.116
eigen(steinquad$asymp.cov, only.values = TRUE)$values 
# largest =  3.297248e+03, smallest = 4.607659e-04

saveRDS(out, file = "GMCBMH_Sitka1988_output_0112001.rds")
