
#### GBR density function ####

# a function to plot the marginal density of the regression coefficient under the GBR prior

g <- function(e, f, beta) {
  function(alpha) {
    alpha * gamma(e + 1/alpha)/gamma(1/alpha) * 1/(abs(beta)^alpha + 2*f)^(e + 1/alpha)
  }
}

dgbr <- function(e1, f1, e2, f2, k1, k2) {
  function(beta) {
    out <- numeric(length(beta))
    for (i in 1:length(beta)) {
      geval1 <- g(e1, f1, beta[i])
      geval2 <- g(e2, f2, beta[i])
      part1 <- (2*f1)^e1/4/gamma(e1) * stats::integrate(geval1, lower = k1, upper = k2)$value
      part2 <- (2*f2)^e2/4/gamma(e2) * stats::integrate(geval2, lower = k1, upper = k2)$value
      out[i] <- part1 + part2
    }
    return(out)
  }
}
