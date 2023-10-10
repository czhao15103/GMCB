\name{b.to.matrix}
\alias{b.to.matrix}
\title{
Reformatting the mean regression coefficients
}
\description{
Converts an \eqn{m x pq} matrix into a list of matrices
}
\usage{
b.to.matrix(b, p, q)
}
\arguments{
  \item{b}{An \eqn{m x pq} matrix, where each row is a vectorized posterior sample of
the mean regression coefficient matrix \emph{B}.}
  \item{p}{scalar, the number of predictors}.
  \item{q}{scalar, the number of responses}.
}
%\details{
%The function \code{gmcb} outputs vec(\emph{B}). For calculating quantities such as 
%}
\value{
A list of length \emph{m} whose elements are \emph{p x q} matrices
}
\author{
\packageAuthor{GMCB}
}

\seealso{
\code{\link{precmat}}
}
\examples{
m <- 1000
p <- 5
q <- 5

b <- matrix(rnorm(m*p*q), nrow = m, ncol = p*q)
b.mat <- b.to.matrix(b.mat, p, q)
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
