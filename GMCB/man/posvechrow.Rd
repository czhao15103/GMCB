\name{pos.vechrow}
\alias{pos.vechrow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Indexing the strictly lower triangular elements of a matrix
}
\description{
A function that returns the index of a strictly lower triangular element of a matrix
in the vector formed by concatenating the strictly lower triangular elements of the
matrix by row.
}
\usage{
pos.vechrow(c, j)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{c}{
row of the element.
}
  \item{j}{
column of the element.
}
}
\details{
Because this function is designed for strictly lower triangular elements,
\code{c} should be at least 2, and \code{j} should be less than \code{c}.
}
\value{
The index corresponding to the specified element in the vector formed by 
concatenating the strictly lower triangular elements of the
matrix by row.
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
\code{\link{vech.row}}
}
\examples{

mat <- matrix(1:16, nrow = 4)

lowertrivec <- vech.row(mat)
index <- pos.vechrow(3, 2)

mat[3,2]
lowertrivec[index]

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
