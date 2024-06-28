\name{init_freqbridge_nocovariates}
\alias{init_freqbridge_nocovariates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initializing the covariance parameters in the absence of covariates with 
frequentist bridge regression
}
\description{
This function assumes that the distribution of the response
has a non-zero mean that does not depend on covariates and 
uses a sequence of frequentist bridge regressions to initialize the covariance 
parameters. 
The regularization parameters and the penalty exponents corresponding to the
mean and covariance are assumed to be the same and selected using k-fold 
cross-validation.
}
\usage{
init_freqbridge_nocovariates(y, alpha.prior, nfolds = 5, nalpha = 10, 
                             nlambda = 20, cores = 4, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response matrix with dimensions \emph{n} by \emph{q}. Must have \eqn{q \ge 2}. 
}
  \item{alpha.prior}{
vector of length 2 giving the minimum and maximum values for the uniform prior
on the penalty exponents.
}
  \item{nfolds}{
folds to use for k-fold cross-validation. 
Default value is 5.
}
  \item{nalpha}{
number of penalty exponent values to consider. 
Default value is 10.
}
  \item{nlambda}{
number of regularization values to consider. 
Default value is 20.
}
  \item{cores}{
number of cores to use for parallelization. 
Set to 1 for no parallelization. Default value is 4.
}
  \item{seed}{
optional scalar passed to \code{set.seed}. Default is \code{NULL}.
}
}
\details{
This function considers \code{nalpha} times \code{nlambda} sets of parameters
for the regularization parameters and penalty exponent parameters in the GMCB
model.  

The candidate penalty exponent values are a sequence of equally spaced
values starting from \code{alpha.prior[1]} with maximal end value 
\code{alpha.prior[2]}. The candidate regularization parameter values are
given by \code{2^((1:nlambda) - 6)}, which was used by Park and Yoon (2011)
in their simulations.
}
\value{
A list including the following components:
  \item{delta.init}{a vector of length \emph{q(q-1)/2}, the initializing values
                    for the covariance parameters with interpretations as 
                    regression coefficients.}
  \item{gamma.init}{a vector of length \emph{q}, the initializing values for 
                    the covariance parameters with interpretations as 
                    regression error variances.}
  \item{gamma.prior}{a vector of length 2, giving the shape and rate 
                     parameters for the inverse Gamma prior on the 
                     covariance parameters with interpretations as 
                     regression error variances.}
  \item{alpha.b.init}{the initializing value for the penalty exponent
                      for the mean.}
  \item{alpha.d.init}{the initializing value for the penalty exponent
                      for the covariance parameters.}
  \item{lambda.init}{a vector of length \emph{q}, the initializing values for the 
                     regularization parameters for the mean.}
  \item{tau.init}{a vector of length \emph{q(q-1)/2}, the initializing values for the 
                  regularization parameters for the covariance parameters.}
  \item{cvresults}{matrix of the k-fold cross-validation mean squared errors
                   for each set of candidate parameter values.}
}
\references{
  Park, C. and Yoon, Y. J. (2011) Bridge regression: Adaptivity and group selection. 
  \emph{Journal of Statistical Planning and Inference} \bold{141} 3506--3519.
  
  Yuzbasi B, Arashi M, Akdeniz F (2020). _rbridge: Restricted Bridge Estimation_. R package version 1.0.2,
  <https://CRAN.R-project.org/package=rbridge>.
}
%\author{
%Bahadir Yuzbasi, Mohammad Arashi, and Fikri Akdeniz wrote the C++ 
%source code for computing the frequentist bridge. Christina Zhao
%wrote the wrapper that adapts the frequentist bridge for 
%initializing the GMCB model.
%}
\note{
It is recommended that \code{nalpha}\eqn{\cdot}\code{nlambda} not be too large
to avoid long computation times. Even with the default values and parallelization, 
runtime can be extensive if \emph{q} is large.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{init_naive}}, \code{\link{nuisance.init}}, \code{\link{ols.init}}
}
\examples{

library(mvtnorm)

n <- 100
q <- 5
rho <- 0.7
y.mats <- matgen(q, rho, type = "AR")
y <- rmvnorm(n, mean = rep(0, q), sigma = y.mats$cov.mat)

alpha.prior <- c(0.5, 4)

initialize <- init_freqbridge_nocovariates(y = y, alpha.prior = alpha.prior,
                                           cores = 2, seed = 1045)

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
