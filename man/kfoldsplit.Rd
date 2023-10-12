\name{kfoldsplit}
\alias{kfoldsplit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
k-fold cross-validation
}
\description{
A function for dividing a set of observations for \emph{k}-fold
cross-validation. The generated folds do not overlap.
}
\usage{
kfoldsplit(n, folds = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of observations.
}
  \item{folds}{
number of folds. Default value is 5.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list of length \code{folds}, with each element containing the indices for that fold.
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

n <- 1000
folds <- kfoldsplit(n)

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
