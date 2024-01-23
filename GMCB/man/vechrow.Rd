\name{vech.row}
\alias{vech.row}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Concatenating the lower triangular elements of a matrix
}
\description{
A function that converts a square matrix to a vector by stacking the
lower-triangular elements by row. 
}
\usage{
vech.row(mat, keep.diag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
matrix, should be square.
}
  \item{keep.diag}{
logical. Default value is \code{FALSE}. If \code{TRUE}, includes the diagonals during concatentation.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
If \code{mat} is a \emph{q} by \emph{q} matrix and \code{keep.diag = FALSE}, a vector of length of \emph{q(q-1)/2}
containing all strictly lower-triangular elements of \code{mat} is returned.
If \code{keep.diag = TRUE}, a vector of length of \emph{q(q+1)/2} is returned.
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
\code{\link{pos.vechrow}}
}
\examples{

mat <- matrix(1:16, nrow = 4)

nodiag <- vech.row(mat)
keepdiag <- vech.row(mat, keep.diag = TRUE)

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
