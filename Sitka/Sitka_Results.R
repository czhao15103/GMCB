
#### Sitka inference ####

library(GMCB)

data <- readRDS("Sitka1988_arranged.rds")

p <- ncol(data$x)
q <- ncol(data$y)

mh0112001 <- readRDS("GMCBMH_Sitka1988_output_0112001.rds")
mh114005 <- readRDS("GMCBMH_Sitka1988_output_114005.rds")
smn0112001 <- readRDS("GMCBSMN_Sitka1988_output_0112001.rds")
smn114005 <- readRDS("GMCBSMN_Sitka1988_output_114005.rds")

mh0112001_b <- b.to.matrix(mh0112001$mcmc$b, p, q)
mh114005_b <- b.to.matrix(mh114005$mcmc$b, p, q)
smn0112001_b <- b.to.matrix(smn0112001$mcmc$b, p, q)
smn114005_b <- b.to.matrix(smn114005$mcmc$b, p, q)

postmean <- list(mh0112001 = Reduce(`+`, mh0112001_b)/length(mh0112001_b),
                 mh114005 = Reduce(`+`, mh114005_b)/length(mh114005_b),
                 smn0112001 = Reduce(`+`, smn0112001_b)/length(smn0112001_b),
                 smn114005 = Reduce(`+`, smn114005_b)/length(smn114005_b))

postmean
postmean$smn0112001

# posterior credible intervals - 95%
postcred <- list(mh0112001 = apply(mh0112001$mcmc$b, 2, quantile, probs = c(0.025, 0.975)),
                 mh114005 = apply(mh114005$mcmc$b, 2, quantile, probs = c(0.025, 0.975)),
                 smn0112001 = apply(smn0112001$mcmc$b, 2, quantile, probs = c(0.025, 0.975)),
                 smn114005 = apply(smn114005$mcmc$b, 2, quantile, probs = c(0.025, 0.975)))
postcred_mat <- lapply(postcred, b.to.matrix, p, q)
for (i in 1:length(postcred_mat)) {
  names(postcred_mat[[i]]) <- c("2.5%", "97.5%")
}
postcred_mat
postcred_mat$smn0112001

# length of the credible intervals
(postcred95length <- lapply(postcred_mat, function(x) x$`97.5%` - x$`2.5%`))
sapply(postcred95length, mean)
#  mh0112001   mh114005 smn0112001  smn114005 
# 0.06314618 0.09237099 0.04916853 0.08918292

mh0112001_omega <- prec.mat(mh0112001$mcmc$delta, mh0112001$mcmc$gamma)
mh114005_omega <- prec.mat(mh114005$mcmc$delta, mh114005$mcmc$gamma)
smn0112001_omega <- prec.mat(smn0112001$mcmc$delta, smn0112001$mcmc$gamma)
smn114005_omega <- prec.mat(smn114005$mcmc$delta, smn114005$mcmc$gamma)

postmean_cov <- list(mh0112001 = Reduce(`+`, mh0112001_omega$cov.mat)/length(mh0112001_omega$cov.mat),
                     mh114005 = Reduce(`+`, mh114005_omega$cov.mat)/length(mh114005_omega$cov.mat),
                     smn0112001 = Reduce(`+`, smn0112001_omega$cov.mat)/length(smn0112001_omega$cov.mat),
                     smn114005 = Reduce(`+`, smn114005_omega$cov.mat)/length(smn114005_omega$cov.mat))
postmean_cov

# smn0112001
blables <- c(bquote(B[11]), bquote(B[12]), bquote(B[13]), bquote(B[14]), bquote(B[15]),
             bquote(B[21]), bquote(B[22]), bquote(B[23]), bquote(B[24]), bquote(B[25]))
deltalab <- c(expression(delta[paste(2, ",", 1, sep = "")]), 
              expression(delta[paste(3, ",", 1, sep = "")]), 
              expression(delta[paste(3, ",", 2, sep = "")]), 
              expression(delta[paste(4, ",", 1, sep = "")]), 
              expression(delta[paste(4, ",", 2, sep = "")]), 
              expression(delta[paste(4, ",", 3, sep = "")]),
              expression(delta[paste(5, ",", 1, sep = "")]), 
              expression(delta[paste(5, ",", 2, sep = "")]), 
              expression(delta[paste(5, ",", 3, sep = "")]),
              expression(delta[paste(5, ",", 4, sep = "")]))

par(mfrow = c(1,2))
for (i in 1:ncol(smn0112001$mcmc$b)) {
  plot(ts(smn0112001$mcmc$b[(nrow(smn0112001$mcmc$b) - 999):nrow(smn0112001$mcmc$b), i]), 
       main = blables[i], ylab = "")
  acf(smn0112001$mcmc$b[,i], lag.max = 100, main = "")
}

par(mfrow = c(1,2))
for (i in 1:ncol(smn0112001$mcmc$delta)) {
  plot(ts(smn0112001$mcmc$delta[(nrow(smn0112001$mcmc$delta) - 999):nrow(smn0112001$mcmc$delta), i]), 
       , main = deltalab[i], ylab = "")
  acf(smn0112001$mcmc$delta[,i], lag.max = 100, main = "")
}

par(mfrow = c(1,2))
for (i in 1:ncol(smn0112001$mcmc$gamma)) {
  plot(ts(smn0112001$mcmc$gamma[(nrow(smn0112001$mcmc$gamma) - 999):nrow(smn0112001$mcmc$gamma), i]), 
       , main = bquote(gamma[.(i)]), ylab = "")
  acf(smn0112001$mcmc$gamma[,i], lag.max = 100, main = "")
}
