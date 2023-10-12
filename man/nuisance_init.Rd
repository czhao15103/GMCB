\name{nuisance.init}
\alias{nuisance.init}
\alias{nuisance.init.meanzero}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initializing the regularization parameters and penalty exponents
}
\description{
User-specified priors are used to generate initializing values for the regularization
parameters and penalty exponents. 
}
\usage{
nuisance.init(p, q, lambda.prior, tau.prior, alpha.prior, 
              b.large = TRUE, delta.large = TRUE)

nuisance.init.meanzero(q, tau.prior, alpha.prior, delta.large = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
number of covariates.
}
  \item{q}{
dimension of response.
}
  \item{lambda.prior}{
either a vector of length 4 or a matrix with dimensions 4 by \code{pq}. 
See `Details.`
}
  \item{tau.prior}{
either a vector of length 4 or a matrix with dimensions 4 by \code{q(q-1)/2}. 
See `Details.`
}
  \item{alpha.prior}{
vector of length 2 giving the minimum and maximum values for the uniform prior
on the penalty exponents.
}
  \item{b.large}{
logical, default value is \code{TRUE}. See `Details.`
}
  \item{delta.large}{
logical, defaults value is \code{TRUE}. See `Details.`
}
}
\details{
The function \code{nuisance.init} assumes a multivariate linear regression model 
for the response, while \code{nuisance.init.meanzero} assumes the mean of the 
response is zero.

The priors on the regularization parameters for the mean regression 
coefficient matrix and for the covariance parameters are two-component Gamma
mixture priors. The mean and variance of the first component are assumed to be 
small, while the mean and variance of the second component are assumed to be large.

For both \code{lambda.prior} and \code{tau.prior}, if the function is provided
with a vector, 
the first two elements provide the shape and rate of the first component of the 
Gamma mixture prior, and the last two elements provide the shape and rate of the
second component of the mixture prior. If the function is provided with a matrix,
each column provides the parameters of the mixture priors for each regularization
parameter, with the first two rows providing the shape and rate of the first
component and the last two rows providing the shape and rate of the second component.

Since large values of the corresponding regression coefficients
are more likely when the regularization parameters are small, initializing values
for the regularization parameters are generated from the first component of 
the mixture when \code{b.large = TRUE} or \code{delta.large = TRUE}. Initializing
values are generated from the second component otherwise.
}
\value{
Both functions require a list that includes the following components:
  \item{lambda.init}{a vector of length \emph{pq}, the initializing values
                     for the regularization parameters associated with the 
                     mean regression coefficient matrix. Only returned if using \code{ols.init}.}
  \item{tau.init}{a vector of length \emph{q(q-1)/2}, the initializing values
                  for the regularization parameters associated with the 
                  covariance parameters.}
  \item{alpha.b.init}{the initializing value for the penalty exponent
                      associated with the mean regression coefficient matrix.
                      Only returned if using \code{ols.init}.}
  \item{alpha.d.init}{the initializing value for the penalty exponent
                      associated with the covariance parameters.}
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
\code{\link{init_freqbridge_nocovariates}}, \code{\link{init_naive}}, 
\code{\link{ols_init}}
}
\examples{

## nuisance.init example
p <- 5
q <- 5

lambda.prior <- c(1, 1, 40, 0.5)
tau.prior <- c(1, 1, 40, 0.5)
alpha.prior <- c(0.5, 4)

initializenuisance <- nuisance.init(p, q, lambda.prior, tau.prior, alpha.prior, 
                                    b.large = TRUE, delta.large = TRUE)
                                    
## nuisance.init.meanzero example
q <- 5

tau.prior <- c(1, 1, 40, 0.5)
alpha.prior <- c(0.5, 4)

initializenuisance <- nuisance.init.meanzero(q, tau.prior, alpha.prior, delta.large = TRUE)

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
