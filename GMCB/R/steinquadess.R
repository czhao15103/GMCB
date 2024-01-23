
# a function for computing the multivariate ESS for the Bayes estimate
# of (B, Omega) under the sum of scalar quadratic loss and Stein's loss

ess_steinquad <- function(b, delta, gamma) {
  m <- nrow(b)
  q <- ncol(gamma)
  p <- ncol(b)/q
  
  stopifnot(ncol(delta) == q*(q-1)/2)
  stopifnot(nrow(delta) == m)
  stopifnot(nrow(gamma) == m)
  
  # create duplication matrix and its Moore-Penrose inverse
  dq <- matrixcalc::duplication.matrix(q)
  dqt.dq <- crossprod(dq)
  dq.mp <- tcrossprod(solve(dqt.dq), dq)
  
  # zero matrices
  zerotopright <- matrix(0, nrow = p*q, ncol = q*(q+1)/2)
  zerobottomleft <- matrix(0, nrow = q*(q+1)/2, ncol = p*q + q*(q + 1)/2)
  
  matrices <- prec.mat(delta, gamma, cov = TRUE)
  omega.mat <- vech.list(matrices$prec.mat) # omega
  omegainv.mat <- vech.list(matrices$cov.mat) # omegainv
  b.matrix <- b.to.matrix(b = b, p = p, q = q)
  bomega <- lapply(1:m, 
                   function(x) b.matrix[[x]] %*% matrices$prec.mat[[x]])
  bomega.mat <- t(sapply(bomega, matrixcalc::vec))
  
  # asymptotic covariance matrix of B\Omega, Omega, and Omega inv
  first.chain <- cbind(bomega.mat, omega.mat, omegainv.mat)
  first.chain.asymptoticcov <- mcmcse::mcse.multi(first.chain, r = 1)$cov
  
  # use delta method with the following matrix function g
  # g(m1, m2, m3) = (m1 %*% m2^{-1}, m3^{-1}) where m1, m2, m3 are pxq, qxq symmetric, and qxq symmetric, respectively
  
  # estimate of m1 = posterior mean of B\Omega
  m1 <- Reduce(`+`, bomega)/m
  
  # estimate of m2^{-1} = inverse of posterior mean of Omega
  m2 <- Reduce(`+`, matrices$prec.mat)/m
  m2inv <- chol2inv(chol(m2))
  
  # estimate of m3^{-1} = Bayes estimate of Omega under Stein's loss
  m3 <- Reduce(`+`, matrices$cov.mat)/m
  m3inv <- chol2inv(chol(m3))
  
  ## constructing the derivative matrix required for delta method
  dvecm1 <- kronecker(m2inv, diag(p))
  dvechm2 <- -kronecker(m2inv, m1 %*% m2inv) %*% dq
  dvechm3 <- -dq.mp %*% kronecker(m3inv, m3inv) %*% dq
  gradg <- rbind(cbind(dvecm1, dvechm2, zerotopright),
                 cbind(zerobottomleft, dvechm3))
  
  # desired asymptotic covariance matrix
  final.asymp.cov <- gradg %*% tcrossprod(first.chain.asymptoticcov, gradg)
  eigs_cov <- eigen(final.asymp.cov, only.values = TRUE)$values # following multiESS in mcmcse (Flegal et al, 2021)
  
  # applying g to the chain, we require sample covariance of (B, \Omega)
  b_omega <- cbind(b, omega.mat)
  b_omega_sampcov <- cov(b_omega)
  
  # compute multivariate ESS - the following code from multiESS in mcmcse (Flegal et al, 2021)
  log.det.var.p <- sum(log(eigen(b_omega_sampcov, symmetric = TRUE, 
                                 only.values = TRUE)$values))
  log.det.covmat.p <- sum(log(eigs_cov))
  ess <- nrow(first.chain) * exp((log.det.var.p - log.det.covmat.p)/ncol(b_omega))
  
  out <- list(multiESS = ess, asymp.cov = final.asymp.cov)
  return(out)
}