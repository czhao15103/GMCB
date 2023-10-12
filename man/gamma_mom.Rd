\name{gamma.mom}
\alias{gamma.mom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Eliciting a prior for gamma
}
\description{
Elicits an inverse Gamma prior using a method of moments estimator.
}
\usage{
gamma.mom(gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gamma}{
vector of initializing values.
}
}
%\details{
%Computes the shape and rate parameter of an inverse Gamma prior
%}
\value{
A vector of length 2. The first component is the shape parameter of 
the prior, and the second is the rate parameter.
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

gamma <- rgamma(20, shape = 3, rate = 2)
gamma.prior <- gamma.mom(gamma)

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
