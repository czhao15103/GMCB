
#### Exponential power distribution ####

library(rmutil)
par(mfrow = c(2, 4))
lambda <- c(1, 5)
alpha <- c(1, 2, 4, 10)
par(mfrow = c(length(lambda), length(alpha)), cex.main = 3)
for (i in 1:length(lambda)) {
  for (j in 1:length(alpha)) {
    f <- alpha[j]/2
    s <- (1/lambda[i])^(2/alpha[j])
    curve(dpowexp(x, m = 0, s = s, f = f), xlim = c(-3, 3), ylim = c(0, 2),
          main = bquote(alpha~"="~.(alpha[j])~","~lambda~"="~.(lambda[i])), ylab = "")
  }
}

#### GBR Prior Plotting ####

# functions to plot the marginal density of the regression coefficient under the GBR prior
g <- function(e, f, beta) {
  function(alpha) {
    alpha * gamma(e + 1/alpha)/gamma(1/alpha) * 1/(abs(beta)^alpha + 2*f)^(e + 1/alpha)
  }
}

f <- function(e1, f1, e2, f2, k1, k2) {
  function(beta) {
    out <- numeric(length(beta))
    for (i in 1:length(beta)) {
      geval1 <- g(e1, f1, beta[i])
      geval2 <- g(e2, f2, beta[i])
      part1 <- (2*f1)^e1/4/gamma(e1) * integrate(geval1, lower = k1, upper = k2)$value
      part2 <- (2*f2)^e2/4/gamma(e2) * integrate(geval2, lower = k1, upper = k2)$value
      out[i] <- part1 + part2
    }
    return(out)
  }
}

#### Plots for section "The GBR Density Function" ####

lambda.prior <- c(0.1, 1, 2, 0.01)
alpha.lb <- c(1e-2, 0.5, 1)
alpha.ub <- c(2, 4, 8)
x <- seq(-4, 4, length.out = 1e4)

# for the typical lower bound of 1
par(mfrow = c(1,1))  
lty <- c(1, 2, 3)
ff <- f(lambda.prior[1], lambda.prior[2], lambda.prior[3], lambda.prior[4],
        alpha.lb[3], alpha.ub[1])
yy <- ff(x)
plot(x, yy, type = "l",xlim = c(-1, 1), ylim = c(0, 5),
     xlab = bquote(beta), ylab = "Density")
for (j in 2:3) {
  ff <- f(lambda.prior[1], lambda.prior[2], lambda.prior[3], lambda.prior[4],
          alpha.lb[3], alpha.ub[j])
  curve(ff(x), type = "l", lty = lty[j], add = TRUE)
}
# legend("topright", legend = c(bquote("k"[2]~"="~.(alpha.ub[1])),
#                               bquote("k"[2]~"="~.(alpha.ub[2])),
#                               bquote("k"[2]~"="~.(alpha.ub[3]))),
#        col = "black",
#        lty = lty)

# for the typical upper bound of 2, value at zero for different lower bounds
alpha.lb <- c(0.01, 0.1, 0.5, 1)
for (i in 1:length(alpha.lb)) {
  ff <- f(lambda.prior[1], lambda.prior[2], lambda.prior[3], lambda.prior[4],
          alpha.lb[i], alpha.ub[1])
  print(ff(0))
}

#### Plots for hyperparameter selection ####

fdense_rec <- f(1, 1, 40, 0.5, 0.5, 4)
fsparse_rec <- f(0.1, 1, 2, 0.01, 0.5, 4)

par(mfrow = c(1,2))
curve(fdense_rec, from = -0.5, to = 0.5, n = 1e4, 
      main = "Comparison of Priors", xlab = bquote(beta),
      ylab = "Density")
curve(fsparse_rec, add = TRUE, n = 1e4, lty = 5)
curve(fdense_rec, from = 5, to = 20, n = 1e4, 
      main = "Tail Behavior", xlab = bquote(beta),
      ylab = "Density")
curve(fsparse_rec, add = TRUE, n = 1e4, lty = 5)

fdense_usual <- f(1, 1, 40, 0.5, 1, 2)
fsparse_usual <- f(0.1, 1, 2, 0.01, 1, 2)

par(mfrow = c(1,2))
curve(fdense_usual, from = -0.5, to = 0.5, n = 1e4, 
      main = "Comparison of Priors", xlab = bquote(beta),
      ylab = "Density")
curve(fsparse_usual, add = TRUE, n = 1e4, lty = 5)
curve(fdense_usual, from = 5, to = 20, n = 1e4, 
      main = "Tail Behavior", xlab = bquote(beta),
      ylab = "Density")
curve(fsparse_usual, add = TRUE, n = 1e4, lty = 5)
