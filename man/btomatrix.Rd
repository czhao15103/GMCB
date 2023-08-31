\name{b.to.matrix}
\alias{b.to.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reformatting the mean regression coefficients
}
\description{
Converts an \eqn{m x pq} matrix into a list of matrices
}
\usage{
b.to.matrix(b, p, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{An \eqn{m x pq} matrix, where each row is a vectorized posterior sample of
the mean regression coefficient matrix \emph{B}.}
  \item{p}{scalar, the number of predictors}.
  \item{q}{scalar, the number of responses}.
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of length \emph{m} whose elements are \emph{p x q} matrices
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
\packageAuthor{GMCB}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
