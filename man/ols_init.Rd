\name{ols.init}
\alias{ols.init}
\alias{ols.covariance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initializing with OLS regression
}
\description{
These functions initialize the mean and covariance parameters using OLS regression.
}
\usage{
ols.init(y, x)

ols.covariance(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response matrix with dimensions \emph{n} by \emph{q}. Must have \eqn{q \ge 2}. 
}
  \item{x}{
covariate matrix with dimensions \emph{n} by \emph{p}.
}
}
\details{
The function \code{ols.init} assumes a multivariate linear regression model for 
the response, while the function \code{ols.covariance} assumes a non-zero mean
that does not depend on covariates. For the OLS estimates to exist, both functions 
require \eqn{q \le n}, and \code{ols.init} also requires \eqn{p < n}.
}
\value{
A list that includes the following components:
  \item{b.init}{a vector of length \emph{pq}, the initializing values
                for the mean regression coefficient matrix. Only included if using \code{ols.init}.}
  \item{delta.init}{a vector of length \emph{q(q-1)/2}, the initializing values
                    for the covariance parameters with interpretations as 
                    regression coefficients.}
  \item{gamma.init}{a vector of length \emph{q}, the initializing values for 
                    the covariance parameters with interpretations as 
                    regression error variances.}
  \item{gamma.prior}{a vector of length two, giving the shape and rate 
                     parameters for the inverse Gamma prior on the 
                     covariance parameters with interpretations as 
                     regression error variances.}
  
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%This function uses the \code{lm}, which was written by the R Core Team,
%to compute the necessary OLS regressions. \packageAuthor{GMCB} 
%wrote the wrapper that adapts OLS regression for
%initializing the GMCB model.
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{init_freqbridge_nocovariates}}, \code{\link{init_naive}}, 
\code{\link{nuisance_init}}
}
\examples{

library(mvtnorm)

## ols.init example
n <- 100
p <- 5
q <- 5

rho <- 0.7

x.mats <- matgen(p, rho, type = "AR")
y.mats <- matgen(q, rho, type = "AR")

x <- rmvnorm(n, mean = rep(0, p), sigma = x.mats$cov.mat)
b <- matrix(rnorm(p*q), nrow = p)
y <- x \%*\% b + rmvnorm(n, mean = rep(0, q), sigma = y.mats$cov.mat)

initialize <- ols.init(y, x)

## ols.covariance example
n <- 100
q <- 5
rho <- 0.7
y.mats <- matgen(q, rho, type = "AR")
y <- rmvnorm(n, mean = rep(0, q), sigma = y.mats$cov.mat)

initialize <- ols.covariance(y)

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
