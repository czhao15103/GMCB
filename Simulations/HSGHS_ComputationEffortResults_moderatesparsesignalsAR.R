
library(mcmcse)
library(GMCB)
library(R.matlab)
library(matrixcalc)

# file location
fileprefix <- "n100p5q5_moderatesparsesignalsAR_ComputationComparison"

p <- 5
q <- 5
m <- 1e5 # number of iterations in each run

clus <- makeCluster(20)
registerDoParallel(clus)

set.seed(112)
ests <- foreach(i = 1:2000, .packages = c("R.matlab", "mcmcse", "matrixcalc", "GMCB")) %dorng% {
	run <- readMat(paste(fileprefix, i, ".mat", sep = ""))
	b.matrix <- lapply(1:m, function(x) run$beta.save[,,x])
	b.vec <- t(sapply(b.matrix, function(x) as.numeric(vec(x))))

	omega <- lapply(1:m, function(x) run$omega.save[,,x])
	sigma <- lapply(omega, solve, tol = 1e-30)
	omega.mat <- vech.list(omega)
	omegainv.mat <- vech.list(sigma)

	#### multivarate ESS for Bayes estimate under L_F
	output.of.interest <- cbind(b.vec, omega.mat)
	ess_postmean <- multiESS(output.of.interest, r = 1) 

	cltcovest.cov <- mcse.multi(output.of.interest, r = 1)$cov
	cltcovest_eigens <- eigen(cltcovest.cov, only.values = TRUE)$values 
	ess_postmean_cltcovest_largesteigen <- cltcovest_eigens[1]
	ess_postmean_cltcovest_smallesteigen <- cltcovest_eigens[length(cltcovest_eigens)]

	#### multivariate ESS for Bayes estimate under L_Q + L_S
	# create duplication matrix and its Moore-Penrose inverse
	dq <- matrixcalc::duplication.matrix(q)
	dqt.dq <- crossprod(dq)
	dq.mp <- tcrossprod(solve(dqt.dq), dq)

	# zero matrices
	zerotopright <- matrix(0, nrow = p*q, ncol = q*(q+1)/2)
	zerobottomleft <- matrix(0, nrow = q*(q+1)/2, ncol = p*q + q*(q + 1)/2)

	bomega <- lapply(1:m, 
                   function(x) b.matrix[[x]] %*% omega[[x]])
	bomega.mat <- t(sapply(bomega, matrixcalc::vec))

	# asymptotic covariance matrix of B\Omega, Omega, and Omega inv
	first.chain <- cbind(bomega.mat, omega.mat, omegainv.mat)
	first.chain.asymptoticcov <- mcmcse::mcse.multi(first.chain, r = 1)$cov

	# use delta method with the following matrix function g
	# g(m1, m2, m3) = (m1 %*% m2^{-1}, m3^{-1}) where m1, m2, m3 are pxq, qxq symmetric, and qxq symmetric, respectively

	# estimate of m1 = posterior mean of B\Omega
	m1 <- Reduce(`+`, bomega)/m

	# estimate of m2^{-1} = inverse of posterior mean of Omega
	m2 <- Reduce(`+`, omega)/m
	m2inv <- chol2inv(chol(m2))

	# estimate of m3^{-1} = Bayes estimate of Omega under Stein's loss
	m3 <- Reduce(`+`, sigma)/m
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
	b_omega_sampcov <- cov(output.of.interest)

	# compute multivariate ESS - the following code from multiESS in mcmcse (Flegal et al, 2021)
	log.det.var.p <- sum(log(eigen(b_omega_sampcov, symmetric = TRUE, 
		                         only.values = TRUE)$values))
	log.det.covmat.p <- sum(log(eigs_cov))
	steinquad_ess <- nrow(first.chain) * exp((log.det.var.p - log.det.covmat.p)/ncol(output.of.interest))

	ess_steinquad_cltcovest_largesteigen <- eigs_cov[1]
	ess_steinquad_cltcovest_smallesteigen <- eigs_cov[length(eigs_cov)]

	return(list(postmean_ess = ess_postmean, postmean_covtr = sum(cltcovest_eigens),
				postmean_largesteigen = ess_postmean_cltcovest_largesteigen,
				postmean_smallesteigen = ess_postmean_cltcovest_smallesteigen,
				steinquad_ess = steinquad_ess, 
				steinquad_covtr = sum(eigs_cov) ,
				steinquad_largesteigen = ess_steinquad_cltcovest_largesteigen,
				steinquad_smallesteigen = ess_steinquad_cltcovest_smallesteigen))
}

stopCluster(clus)

saveRDS(ests, file = "ComputationEffortResults_moderatesparsesignalsAR.rds")
