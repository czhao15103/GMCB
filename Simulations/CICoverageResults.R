
#### n = 100, p = 5, q = 5, CI Coverage ####

library(GMCB)

data_dp <- readRDS("Data_n100p5q5_DP2002_3A_XJ2019_1.rds")
data_mod <- readRDS("Data_n100p5q5_moderatesparsesignalsAR.rds")
data_small <- readRDS("Data_n100p5q5_smallsignalsCS.rds")

p <- 5
q <- 5

# true covariance/precision matrix for mod and small
rho <- 0.7
mats_mod <- matgen(q, rho, type = "AR")
mats_small <- matgen(q, rho, type = "CS")

# vec/vech of true parameters
dp.b.vec <- matrixcalc::vec(data_dp$b)
dp.omega.vech <- matrixcalc::vech(data_dp$y.prec)
mod.b.vec <- matrixcalc::vec(data_mod$b)
mod.omega.vech <- matrixcalc::vech(mats_mod$prec.mat)
small.b.vec <- matrixcalc::vec(data_small$b)
small.omega.vech <- matrixcalc::vech(mats_small$prec.mat)

mhcoverage_dp_b <- matrix(nrow = 2000, ncol = p*q)
mhcoverage_dp_omega <- matrix(nrow = 2000, ncol = q*(q+1)/2)
mhcoverage_mod_b <- matrix(nrow = 2000, ncol = p*q)
mhcoverage_mod_omega <- matrix(nrow = 2000, ncol = q*(q+1)/2)
mhcoverage_small_b <- matrix(nrow = 2000, ncol = p*q)
mhcoverage_small_omega <- matrix(nrow = 2000, ncol = q*(q+1)/2)

smncoverage_dp_b <- matrix(nrow = 2000, ncol = p*q)
smncoverage_dp_omega <- matrix(nrow = 2000, ncol = q*(q+1)/2)
smncoverage_mod_b <- matrix(nrow = 2000, ncol = p*q)
smncoverage_mod_omega <- matrix(nrow = 2000, ncol = q*(q+1)/2)
smncoverage_small_b <- matrix(nrow = 2000, ncol = p*q)
smncoverage_small_omega <- matrix(nrow = 2000, ncol = q*(q+1)/2)

start_mh_dp <- "GMCBMH_DP2002_3A_XJ2019_1_CICoverage"
start_smn_dp <- "GMCBSMN_DP2002_3A_XJ2019_1_CICoverage"

start_mh_mod <- "GMCBMH_moderatesparseAR_CICoverage"
start_smn_mod <- "GMCBSMN_moderatesparseAR_CICoverage"

start_mh_small <- "GMCBMH_smallsignalsCS_CICoverage"
start_smn_small <- "GMCBSMN_smallsignalsCS_CICoverage"

for (i in 1:2000) {
	cat("i = ", i, "\n")
	
	mh_cred_dp <- readRDS(paste(start_mh_dp, i, ".rds", sep = ""))
	mh_cred_mod <- readRDS(paste(start_mh_mod, i, ".rds", sep = ""))
	mh_cred_small <- readRDS(paste(start_mh_small, i, ".rds", sep = ""))

	smn_cred_dp <- readRDS(paste(start_smn_dp, i, ".rds", sep = ""))
	smn_cred_mod <- readRDS(paste(start_smn_mod, i, ".rds", sep = ""))
	smn_cred_small <- readRDS(paste(start_smn_small, i, ".rds", sep = ""))

	for (j in 1:(p*q)) {
		mhcoverage_dp_b[i,j] <- dp.b.vec[j] >= mh_cred_dp$b.credible[1,j] & dp.b.vec[j] <= mh_cred_dp$b.credible[2,j]

		mhcoverage_mod_b[i,j] <- mod.b.vec[j] >= mh_cred_mod$b.credible[1,j] & mod.b.vec[j] <= mh_cred_mod$b.credible[2,j]

		mhcoverage_small_b[i,j] <- small.b.vec[j] >= mh_cred_small$b.credible[1,j] & small.b.vec[j] <= mh_cred_small$b.credible[2,j]

		smncoverage_dp_b[i,j] <- dp.b.vec[j] >= smn_cred_dp$b.credible[1,j] & dp.b.vec[j] <= smn_cred_dp$b.credible[2,j]

		smncoverage_mod_b[i,j] <- mod.b.vec[j] >= smn_cred_mod$b.credible[1,j] & mod.b.vec[j] <= smn_cred_mod$b.credible[2,j]

		smncoverage_small_b[i,j] <- small.b.vec[j] >= smn_cred_small$b.credible[1,j] & small.b.vec[j] <= smn_cred_small$b.credible[2,j]
	}

	for (j in 1:(q*(q+1)/2)) {
		mhcoverage_dp_omega[i,j] <- dp.omega.vech[j] >= mh_cred_dp$omega.credible[1,j] & dp.omega.vech[j] <= mh_cred_dp$omega.credible[2,j]

		mhcoverage_mod_omega[i,j] <- mod.omega.vech[j] >= mh_cred_mod$omega.credible[1,j] & mod.omega.vech[j] <= mh_cred_mod$omega.credible[2,j]

		mhcoverage_small_omega[i,j] <- small.omega.vech[j] >= mh_cred_small$omega.credible[1,j] & small.omega.vech[j] <= mh_cred_small$omega.credible[2,j]

		smncoverage_dp_omega[i,j] <- dp.omega.vech[j] >= smn_cred_dp$omega.credible[1,j] & dp.omega.vech[j] <= smn_cred_dp$omega.credible[2,j]

		smncoverage_mod_omega[i,j] <- mod.omega.vech[j] >= smn_cred_mod$omega.credible[1,j] & mod.omega.vech[j] <= smn_cred_mod$omega.credible[2,j]

		smncoverage_small_omega[i,j] <- small.omega.vech[j] >= smn_cred_small$omega.credible[1,j] & small.omega.vech[j] <= smn_cred_small$omega.credible[2,j]
	}

}

b.coverage <- data.frame(mh_dp_b = colMeans(mhcoverage_dp_b), 
		   mh_mod_b = colMeans(mhcoverage_mod_b),
		   mh_small_b = colMeans(mhcoverage_small_b), 
		   smn_dp_b = colMeans(smncoverage_dp_b),
		   smn_mod_b = colMeans(smncoverage_mod_b), 
		   smn_small_b = colMeans(smncoverage_small_b))

omega.coverage <- data.frame(mh_dp_omega = colMeans(mhcoverage_dp_omega),
		   mh_mod_omega = colMeans(mhcoverage_mod_omega),
		   mh_small_omega = colMeans(mhcoverage_small_omega),
		   smn_dp_omega = colMeans(smncoverage_dp_omega),
		   smn_mod_omega = colMeans(smncoverage_mod_omega),
		   smn_small_omega = colMeans(smncoverage_small_omega))

colMeans(b.coverage)

colMeans(omega.coverage)

# sparse b 
cbind(dp.b.vec, colMeans(mhcoverage_dp_b), colMeans(smncoverage_dp_b))

cbind(mod.b.vec, colMeans(mhcoverage_mod_b), colMeans(smncoverage_mod_b))


# dense b
cbind(small.b.vec, colMeans(mhcoverage_small_b), colMeans(smncoverage_small_b))

# sparse omega
cbind(mod.omega.vech, colMeans(mhcoverage_mod_omega), colMeans(smncoverage_mod_omega))
# coverage quite good for true signals in sparse omega

# dense omega
cbind(dp.omega.vech, colMeans(mhcoverage_dp_omega), colMeans(smncoverage_dp_omega))
cbind(small.omega.vech, colMeans(mhcoverage_small_omega), colMeans(smncoverage_small_omega))

