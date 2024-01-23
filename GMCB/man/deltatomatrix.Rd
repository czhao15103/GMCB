\name{delta.to.matrix}
\alias{delta.to.matrix}
\title{
Reformatting the distinct elements of a lower triangular matrix with ones on the
diagonal into a matrix
}
\description{
This function converts the output from \code{gmcb} for the covariance parameters
\eqn{\delta} from a matrix into a list of matrices
}
\usage{
delta.to.matrix(delta, q)
}
\arguments{
  \item{delta}{matrix with \code{q(q-1)/2} columns. Rows are assumed to correspond to individual samples.}
  \item{q}{the number of responses.}
}
\details{
The modified Cholesky decomposition of a \eqn{q \times q} precision matrix \eqn{\Omega} 
is given by
\deqn{
\Omega = T'D^{-1}T
}
where \emph{T} denotes a lower triangular matrix with ones on the diagonal and
\emph{D} denotes a diagonal matrix with positive diagonal elements. 

Denote the strictly lower triangular elements of the \emph{c}th row of \emph{T}
by \eqn{-\delta_c}, which has length \eqn{c-1}, for \eqn{c = 2,\ldots,q}. Define
\deqn{\delta = (\delta_2,\delta_3',\ldots,\delta_q')'} 
The function \code{gmcb} 
returns posterior samples of \eqn{\delta} in a matrix with \eqn{q(q-1)/2} columns.
However, computation of Bayes estimators of \eqn{\Omega} require rearrangement
of these posterior samples into matrices.
}
\value{
A list of \code{q(q-1)/2} lower triangular
matrices with ones on the diagonals.
}
%\author{
%Christina Zhao
%}

\seealso{
\code{\link{prec.mat}}
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

t.list <- delta.to.matrix(out$mcmc$delta, q)

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
