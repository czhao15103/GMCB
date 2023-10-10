\name{ess_steinquad}
\alias{ess_steinquad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate effective sample size for a certain Bayes estimator
}
\description{
Calculates the multivariate effective sample size for the mean regression coefficient matrix and the error precision matrix under scalar quadratic loss and Stein's loss
}
\usage{
ess_steinquad(b, delta, gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{the matrix of vectorized samples of the mean regression coefficient matrix}
  \item{delta}{the matrix of the samples of the distinct covariance regression parameters}
  \item{gamma}{the matrix of the samples of the }
}
\details{
The scalar quadratic loss (Yuasa & Kubokawa, 2021) for an estimate \eqn{\tilde{B}} is given by
\deqn{
  L_S(\tilde{B},B) = tr((\tilde{B}-B)\Omega(\tilde{B}-B)')
}
where \eqn{B'} denotes the transpose of \eqn{B}. Stein's loss for an estimate \eqn{\tilde{\Omega}} is given by
\deqn{
  L_Q(\tilde{\Omega}, \Omega) = tr(\tilde{\Omega}\Omega^{-1})-\log|\tilde{\Omega}\Omega^{-1}|-q
}
where \eqn{q} is the dimension of \eqn{\Omega}.

The Bayes estimators for \eqn{(B,\Omega)} under these loss functions (Yang & Berger, 1994; Yuasa & Kubokawa, 2021)  are
\deqn{
  (\hat{B}_Q,\hat{\Omega}_S) =(E(B\Omega|Y,X)[E(\Omega|Y,X)]^{-1},\  [E(\Omega^{-1}|Y,X)]^{-1}) 
}
This function uses the delta method to compute the multivariate effective sample 
size based on the method of Vats et al. (2019).
}
\value{
A list with the following elements:
  \item{multiESS}{the multivariate effective sample size}
  \item{asymp.cov}{the estimated asymptotic covariance matrix}
}
\references{
Flegal, J. M., Hughes, J., Vats, D., & Dai, N. (2020). Mcmcse: Monte carlo standard
errors for mcmc [R package version 1.4-1]. Riverside, CA, Denver, CO, Coventry,
UK, Minneapolis, MN.

Vats, D., Flegal, J. M., & Jones, G. L. (2019). Multivariate output analysis for Markov
chain Monte Carlo. \emph{Biometrika}, 106 (2), 321–337. \url{https://doi.org/10.1093/biomet/asz002}

Yang, R., & Berger, J. O. (1994). Estimation of a covariance matrix using the reference
prior. \emph{The Annals of Statistics}, 22 (3), 1195–1211.

Yuasa, R., & Kubokawa, T. (2021). Generalized Bayes estimators with closed forms for
the normal mean and covariance matrices.
}
\author{
\packageAuthor{GMCB}
}

\examples{

library(mvtnorm)

n <- 100
p <- 5
q <- 5

rho <- 0.7

x.mats <- matgen(p, rho, type = "AR")
y.mats <- matgen(q, rho, type = "AR")

x <- rmvnorm(n, mean = rep(0, p), sigma = x.mats$cov.mat)
b <- matrix(rnorm(p*q), nrow = p)
y <- x %*% b + rmvnorm(n, mean = rep(0, q), sigma = y.mats$cov.mat)

out <- gmcb(y, x)

ess <- ess_steinquad(out$mcmc$b, out$mcmc$delta, out$mcmc$gamma)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
