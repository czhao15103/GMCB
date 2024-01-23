
#### dense B, CS Sigma empirical risk, p = 5, q = 5, n = 100 ####

library(GMCB)

# settings used
data <- readRDS("Data_n100p5q5_smallsignalsCS.rds")
y <- data$y # non-centered
x <- data$x # non-standardized
b <- data$b # true b matrix

# simulation settings
p <- ncol(x)
q <- ncol(y[[1]])
n <- nrow(y[[1]])
rho <- 0.7

# number of replications
iter <- length(y)

# generate true covariance/precision matrix
mats <- matgen(q, rho, type = "CS")

mh <- readRDS("gmcb_mh_n100p5q5_smallsignalsCS.rds")
mSSL <- readRDS("mSSL_n100p5q5_smallsignalsCS.rds")
smn <- readRDS("gmcb_smn_n100p5q5_smallsignalsCS.rds")

mh$time 
smn$time 
mSSL$time 

#### Frobenius loss ####
hsghs <- read.csv("hsghs_n100p5q5_smallsignalsCS.csv")

# Frobenius loss for B
loss.postmeanb.mh <-  sapply(1:iter,
                             function(x) sum((mh$ests[[x]]$b.postmean - b)^2))

loss.postmeanb.smn <-  sapply(1:iter, 
                             function(x) sum((smn$ests[[x]]$b.postmean - b)^2))

loss.b.mSSL_DCPE <-  sapply(1:iter, 
                             function(x) sum((mSSL$ests[[x]]$b.dcpe - b)^2))

loss.b.mSSL_DPE <-  sapply(1:iter, 
                             function(x) sum((mSSL$ests[[x]]$b.dpe - b)^2))

# Frobenius loss for Omega
loss.omega.mSSL_DCPE <-  sapply(1:iter, 
                                function(x) sum((mSSL$ests[[x]]$omega.dcpe - mats$prec.mat)^2))

loss.omega.mSSL_DPE <-  sapply(1:iter, 
                               function(x) sum((mSSL$ests[[x]]$omega.dpe - mats$prec.mat)^2))

loss.postmeanomega.mh <-  sapply(1:iter,
                             function(x) sum((mh$ests[[x]]$omega.postmean - mats$prec.mat)^2))

loss.postmeanomega.smn <-  sapply(1:iter, 
                              function(x) sum((smn$ests[[x]]$omega.postmean - mats$prec.mat)^2))

loss.steinquadb.mh <-  sapply(1:iter,
                             function(x) sum((mh$ests[[x]]$b.quad - b)^2))

loss.steinquadb.smn <-  sapply(1:iter, 
                              function(x) sum((smn$ests[[x]]$b.quad - b)^2))


loss.steinquadomega.mh <-  sapply(1:iter,
                             function(x) sum((mh$ests[[x]]$omega.stein - mats$prec.mat)^2))

loss.steinquadomega.smn <-  sapply(1:iter, 
                                  function(x) sum((smn$ests[[x]]$omega.stein - mats$prec.mat)^2))

# posterior mean loss contribution
postmean <- data.frame(mh.b = (loss.postmeanb.mh), mh.omega = (loss.postmeanomega.mh),
                       smn.b = (loss.postmeanb.smn), smn.omega = (loss.postmeanomega.smn),
                       mSSL_DCPE.b = (loss.b.mSSL_DCPE), mSSL_DCPE.omega = (loss.omega.mSSL_DCPE),
                       mSSL_DPE.b = (loss.b.mSSL_DPE), mSSL_DPE.omega = (loss.omega.mSSL_DPE),
                       hsghs.b = (hsghs$f_losspostmeanb), hsghs.omega = (hsghs$f_loss_postmeanomega))

steinquad <- data.frame(mh.b = (loss.steinquadb.mh), mh.omega = (loss.steinquadomega.mh),
                        smn.b = (loss.steinquadb.smn), smn.omega = (loss.steinquadomega.smn),
                        mSSL_DCPE.b = (loss.b.mSSL_DCPE), mSSL_DCPE.omega = (loss.omega.mSSL_DCPE),
                        mSSL_DPE.b = (loss.b.mSSL_DPE), mSSL_DPE.omega = (loss.omega.mSSL_DPE),
                        hsghs.b = (hsghs$f_loss_quadb), hsghs.omega = (hsghs$f_loss_steinomega))

(loss_f <- apply(postmean, 2, mean))
(se_f <- apply(postmean, 2, function(x) sd(x)/sqrt(iter)))

(loss_sq <- apply(steinquad, 2, mean))
(se_sq <- apply(steinquad, 2, function(x) sd(x)/sqrt(iter)))

#### Scalar quadratic and Stein's loss ####

library(R.matlab)
hsghs.sq <- readMat("n100p5q5_smallsignalsCS.mat")

scalarquadloss <- function(b.est, b.actual, omega.actual) {
  b.diff <- b.est - b.actual
  return(sum(diag(b.diff %*% omega.actual %*% t(b.diff))))
}

kl.loss <- function(omega.est, cov.actual) {
  q <- ncol(omega.est)
  mat.prod <- omega.est %*% cov.actual
  return(sum(diag(mat.prod)) - as.numeric(determinant(mat.prod, logarithm = TRUE)$modulus) - q)
}

# Scalar quadratic loss for B
sqloss.postmeanb.mh <-  sapply(1:iter,
                             function(x) scalarquadloss(mh$ests[[x]]$b.postmean, b, mats$prec.mat))

sqloss.postmeanb.smn <-  sapply(1:iter, 
                              function(x) scalarquadloss(smn$ests[[x]]$b.postmean, b, mats$prec.mat))

sqloss.b.mSSL_DCPE <-  sapply(1:iter, 
                            function(x) scalarquadloss(mSSL$ests[[x]]$b.dcpe, b, mats$prec.mat))

sqloss.b.mSSL_DPE <-  sapply(1:iter, 
                           function(x) scalarquadloss(mSSL$ests[[x]]$b.dpe, b, mats$prec.mat))

sqloss.postmeanb.hsghs <- sapply(1:iter,
                       function(x) scalarquadloss(hsghs.sq$beta.mean[x][[1]][[1]], b, mats$prec.mat))

sqloss.steinquadb.mh <-  sapply(1:iter,
                              function(x) scalarquadloss(mh$ests[[x]]$b.quad, b, mats$prec.mat))

sqloss.steinquadb.smn <-  sapply(1:iter, 
                               function(x) scalarquadloss(smn$ests[[x]]$b.quad, b, mats$prec.mat))

sqloss.steinquadb.hsghs <- sapply(1:iter,
                                  function(x) scalarquadloss(hsghs.sq$beta.quad[x][[1]][[1]], b, mats$prec.mat))

# Stein's loss for Omega
sloss.omega.mSSL_DCPE <-  sapply(1:iter, 
                                function(x) kl.loss(mSSL$ests[[x]]$omega.dcpe, mats$cov.mat))

sloss.omega.mSSL_DPE <-  sapply(1:iter, 
                               function(x) kl.loss(mSSL$ests[[x]]$omega.dpe, mats$cov.mat))

sloss.postmeanomega.mh <-  sapply(1:iter,
                                 function(x) kl.loss(mh$ests[[x]]$omega.postmean, mats$cov.mat))

sloss.postmeanomega.smn <-  sapply(1:iter, 
                                  function(x) kl.loss(smn$ests[[x]]$omega.postmean, mats$cov.mat))

sloss.postmeanomega.hsghs <-  sapply(1:iter, 
                                   function(x) kl.loss(hsghs.sq$omega.mean[x][[1]][[1]], mats$cov.mat))

sloss.steinquadomega.mh <-  sapply(1:iter,
                                  function(x) kl.loss(mh$ests[[x]]$omega.stein, mats$cov.mat))

sloss.steinquadomega.smn <-  sapply(1:iter, 
                                   function(x) kl.loss(smn$ests[[x]]$omega.stein, mats$cov.mat))

sloss.steinquadomega.hsghs <-  sapply(1:iter, 
                                     function(x) kl.loss(hsghs.sq$omega.stein[x][[1]][[1]], mats$cov.mat))

# posterior mean loss contribution
postmean.sq <- data.frame(mh.b = (sqloss.postmeanb.mh), mh.omega = (sloss.postmeanomega.mh),
                       smn.b = (sqloss.postmeanb.smn), smn.omega = (sloss.postmeanomega.smn),
                       mSSL_DCPE.b = (sqloss.b.mSSL_DCPE), mSSL_DCPE.omega = (sloss.omega.mSSL_DCPE),
                       mSSL_DPE.b = (sqloss.b.mSSL_DPE), mSSL_DPE.omega = (sloss.omega.mSSL_DPE),
                       hsghs.b = sqloss.postmeanb.hsghs, hsghs.omega = sloss.postmeanomega.hsghs)

steinquad.sq <- data.frame(mh.b = (sqloss.steinquadb.mh), mh.omega = (sloss.steinquadomega.mh),
                        smn.b = (sqloss.steinquadb.smn), smn.omega = (sloss.steinquadomega.smn),
                        mSSL_DCPE.b = (sqloss.b.mSSL_DCPE), mSSL_DCPE.omega = (sloss.omega.mSSL_DCPE),
                        mSSL_DPE.b = (sqloss.b.mSSL_DPE), mSSL_DPE.omega = (sloss.omega.mSSL_DPE),
                        hsghs.b = sqloss.steinquadb.hsghs, hsghs.omega = sloss.steinquadomega.hsghs)

(loss_f <- apply(postmean.sq, 2, mean))
(se_f <- apply(postmean.sq, 2, function(x) sd(x)/sqrt(iter)))

(loss_sq <- apply(steinquad.sq, 2, mean))
(se_sq <- apply(steinquad.sq, 2, function(x) sd(x)/sqrt(iter)))

