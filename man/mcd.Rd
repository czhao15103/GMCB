\name{mcd}
\alias{mcd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing the modified Cholesky decomposition
}
\description{
Computes the modified Cholesky decomposition of the inverse of the input
covariance matrix.
}
\usage{
mcd(cov.mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.mat}{
a covariance matrix.
}
}
\details{
Suppose that \eqn{\Sigma} is a covariance matrix, and let \emph{R} denote its
lower triangular Cholesky factor. Let \emph{D} be the diagonal matrix with
diagonal elements the squared diagonal elements of \emph{R}. Define
\deqn{
L = RD^{-1/2}
}
Then
\deqn{
\Sigma = LDL'
}
Define \eqn{T = L^{-1}}. Then the modified Cholesky decomposition of \eqn{\Sigma^{-1}}
is 
\deqn{
\Sigma^{-1} = T'D^{-1}T
}
}
\value{
A list including the following components:
  \item{t}{the matrix \emph{T}, a lower triangular matrix with ones on the diagonal.}
  \item{d}{the matrix \emph{D}, a diagonal matrix with postive diagonal elements.}
}
\references{
  Pourahmadi, M. (2011) Covariance Estimation: The GLM and Regularization 
  Perspectives. \emph{Statistical Science} \bold{26}:3 369--387.
}
%\author{
%Christina Zhao
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

mats <- matgen(p = 10, rho = 0.7, type = "CS")
mcd.cs <- mcd(mats$cov.mat)

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
