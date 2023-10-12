\name{b.to.matrix}
\alias{b.to.matrix}
\title{
Reformatting the mean regression coefficients
}
\description{
This function converts the output from \code{gmcb} for the mean regression coefficient matrix
\emph{B} from a matrix into a list of matrices
}
\usage{
b.to.matrix(b, p, q)
}
\arguments{
  \item{b}{matrix with \code{pq} columns. Rows are assumed to correspond to individual samples.}
  \item{p}{number of predictors.}
  \item{q}{number of responses.}
}
\details{
When covariates are provided, the function \code{gmcb} assumes the following
distribution of the \emph{q}-variate response \emph{Y} with \emph{p}-variate
predictor \emph{X}:
\deqn{
Y \sim N(B'X,\Omega^{-1})
}
where \emph{B} is a \eqn{p \times q} matrix, and \eqn{\Omega} 
is a positive-definite \eqn{q \times q} matrix. Posterior samples of \emph{B} 
are vectorized and saved as the rows of a matrix with \emph{pq} columns.
}
\value{
A list of \code{p} by \code{q} matrices
}
%\author{
%Christina Zhao
%}

\seealso{
\code{\link{precmat}}
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

b.list <- b.to.matrix(out$mcmc$b, p, q)

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
