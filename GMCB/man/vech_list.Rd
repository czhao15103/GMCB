\name{vech.list}
\alias{vech.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtaining the vech of a list of matrices
}
\description{
A function that turns a list of length \emph{m} with elements \emph{q} by \emph{q}
matrices into a matrix with dimensions \emph{m} by \emph{q(q+1)/2}, where each row
of the matrix is the vech of the corresponding matrix in the list.
}
\usage{
vech.list(list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list}{
a list of length \emph{m} with elements \emph{q} by \emph{q}
matrices
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A matrix of dimension \emph{m} by \emph{q(q+1)/2}.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
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

preconly <- prec.mat(out$mcmc$delta, out$mcmc$gamma, cov = FALSE)

vechsamples <- vech.list(preconly)

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
