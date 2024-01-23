\name{prec.mat}
\alias{prec.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reformatting the output into posterior samples of the covariance/precision matrix
}
\description{
A function that turns the parameters of the modified Cholesky decomposition
into the \emph{q} by \emph{q} covariance/precision matrix.
}
\usage{
prec.mat(delta, gamma, cov = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
matrix with \emph{q(q-1)/2} columns.
}
  \item{gamma}{
matrix with \emph{q} columns.
}
  \item{cov}{
logical specifying whether to return posterior samples of the covariance matrix. Default value is \code{TRUE}.
}
}
\details{
The modified Cholesky decomposition of a \emph{q} by \emph{q} precision matrix \eqn{\Omega} 
is given by
\deqn{
\Omega = T'D^{-1}T
}
where \emph{T} denotes a lower triangular matrix with ones on the diagonal and
\emph{D} denotes a diagonal matrix with positive diagonal elements. 
Each row of \code{delta} should contain a posterior sample of the strictly
lower-triangular elements of \emph{T}, concatenated row by row. Each row of
\code{gamma} should contain a posterior sample of the diagonal elements of \emph{D}.
}
\value{
Either a list of posterior samples for the precision matrix or a list of two lists, 
one containing samples of the precision matrix, the other the samples of the covariance matrix.
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%\packageAuthor{GMCB}
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gmcb}}
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
y <- x \%*\% b + rmvnorm(n, mean = rep(0, q), sigma = y.mats$cov.mat)

priors <- list(lambda = c(1, 1, 40, 0.5), tau = c(1, 1, 40, 0.5),
               alpha = c(0.5, 4))
rw.scale <- list(b = 0.1, delta = 0.1, alpha.b = 0.5, alpha.d = 0.5)

out <- gmcb(y = y, x = x, priors = priors, rw.scale = rw.scale)

## posterior samples of the precision matrix only
preconly <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = FALSE)

## posterior samples of the precision matrix and covariance matrix
preccov <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = TRUE)

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
