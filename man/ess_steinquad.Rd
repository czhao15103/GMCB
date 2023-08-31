\name{ess_steinquad}
\alias{ess_steinquad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate effective sample size for the Bayes estimate under scalar quadratic loss and Stein's loss
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ess_steinquad(b, delta, gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{the matrix of vectorized samples of the mean regression coefficient matrix}
  \item{delta}{the matrix of the samples of the distinct covariance regression parameters}
  \item{gamma}{the matrix of the samples of the }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the following elements:
  \item{multiESS}{the multivariate effective sample size}
  \item{asymp.cov}{the estimated asymptotic covariance matrix}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
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
