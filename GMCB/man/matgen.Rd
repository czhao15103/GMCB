\name{matgen}
\alias{matgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating compound symmetric, AR(1), or diagonal matrices
}
\description{
A function for quickly generating compound symmetric, AR(1), and diagonal matrices,
as well as the corresponding inverse matrix.
}
\usage{
matgen(p, rho, sigma = 1, type = c("AR", "CS", "cI"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
dimension of the matrix to be generated.
}
  \item{rho}{
correlation when \code{type = AR} or \code{type = CS}. Ignored when
\code{type = cI}.
}
  \item{sigma}{
common variance if \code{type = AR} or \code{type = CS}, 
the diagonal value if \code{type = cI}. Default value is \code{sigma = 1}.
}
  \item{type}{
type of matrix to generate, one of \code{AR}, \code{CS}, or \code{cI}.
Defaults to \code{AR}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list including the following components:
  \item{cov.mat}{\code{p} by \code{p} matrix of the specified type.}
  \item{prec.mat}{\code{p} by \code{p} matrix, inverse of the specified covariance matrix.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
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

p <- 10
rho <- 0.7

mats <- matgen(p, rho, type = "AR")
mats <- matgen(p, rho, type = "CS")
mats <- matgen(p = p, sigma = 3, type = "cI")

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
