\name{gmcb}
\alias{gmcb}
\alias{gmcb.matrix}
\alias{gmcb.gbridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simultaneous mean and covariance estimation under the generalized mean-covariance bridge prior
}
\description{
Runs either the GMCB-MH or GMCB-SMN algorithm on input data and returns the posterior samples.
}
\usage{
gmcb(y, x = NULL, initial.values = NULL, priors = NULL, rw.scale, iter = 1e4, 
     algorithm = c("MH", "SMN"), debug = FALSE)
                 
\method{gmcb}{matrix}(y, x = NULL, initial.values = NULL, priors = NULL, rw.scale, 
                      iter = 1e4, algorithm = c("MH", "SMN"), debug = FALSE)

\method{gmcb}{gbridge}(y, x = NULL, initial.values = NULL, priors = NULL, rw.scale, 
                       iter = 1e4, algorithm = c("MH", "SMN"), debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{
response matrix with dimensions \eqn{n \times q}. Must have \eqn{q \ge 2}. 
}
  \item{x}{
optional covariate matrix with dimensions \eqn{n \times p}. If not provided, function will use an intercept of ones.
}
  \item{initial.values}{
optional list with named elements. If \code{algorithm = MH}, must include named elements
\code{b}, \code{delta}, \code{gamma}, \code{alpha.b}, \code{alpha.d}. If \code{algorithm = SMN}, must include named elements \code{b}, \code{delta}, \code{gamma}, \code{lambda}, \code{tau}, \code{alpha.b}, \code{alpha.d}.
}
  \item{priors}{
list with named elements \code{lambda}, \code{tau}, \code{alpha}. May also include
optional named element \code{gamma}.
}
  \item{rw.scale}{
list with named elements giving the standard deviations of the normal proposals used for the Metropolis-Hastings updates. If \code{algorithm = MH}, must include named elements
\code{b}, \code{delta}, \code{alpha.b}, \code{alpha.d}. If \code{algorithm = SMN}, must include named elements \code{alpha.b}, \code{alpha.d}
}
  \item{iter}{
number of iterations. Default value is \code{1e4}.
}
  \item{debug}{
logical. Default value is \code{FALSE}. If \code{TRUE}, returns the random variates
generated during sampling as described below. 
}
  \item{algorithm}{
specifies the algorithm to run. Either \code{MH} or \code{SMN}.
}
}
\details{
  Let \emph{Y} be a vector of length \emph{q}, \emph{X} a vector length \emph{p}, \emph{B} a \eqn{p \times q}
  matrix, and \eqn{\Omega} a positive-definite \eqn{q \times q} matrix. Suppose that the distribution of
  \emph{Y} is given by \eqn{N(B' X,\Omega^{-1})}, where \emph{B'} denotes the transpose of \emph{B}. 
  
  Define the modified Cholesky decomposition of \eqn{\Omega^{-1}} to be
  \deqn{
    \Omega = T'D^{-1}T
  }
  where \emph{T} is a unit lower triangular matrix and \emph{D} is a diagonal
  matrix with positive entries. By reparameterizing the likelihood of \emph{Y} using
  this modified Cholesky decomposition, the problem of estimating \eqn{(B,\Omega)}
  can be turned into estimating the regression coefficients and error variances of a
  sequence of \emph{q} regressions.
  
  The generalized mean-covariance bridge (GMCB) prior uses the generalized bridge prior proposed
  by Xiang (2020) to model the parameters of the \emph{q} regressions. The posterior distribution
  under the GMCB prior can be sampled using either the GMCB-MH algorithm or the GMCB-SMN algorithm.
  
  The details of the hierarchical model and the algorithms are included in the package 
  vignette \emph{GMCB}, which can be accessed with 
  \code{vignette("GMCB", "GMCB")}.

  This function employs several ideas discussed in Geyer (2011). Users can start runs that pick up
  where previous runs ended by passing an object of class \code{gbridge} to the function, and the 
  \code{debug} output is designed to allow "white box" testing as described in section 1.16 of 
  Geyer (2011).
}

\value{
An object of class \code{gbridge}, which is a list including the following components:
  \item{mcmc}{a list with named elements \code{b}, \code{delta}, \code{gamma},
              \code{lambda}, \code{tau}, \code{alpha_b}, and \code{alpha_d}.}
  \item{acceptances}{a list giving the acceptance rates for the Metropolis-Hastings
                     updates. Has length 4 if \code{algorithm = "MH"} and length 2
                     if \code{algorithm = "SMN"}.}
  \item{start.seed}{the value of \code{.Random.seed} before sampling began.}
  \item{end.seed}{the value of \code{.Random.seed} after sampling ended.}
  \item{timing}{a character vector of length 2, giving the time used for 
                pre-sampling computations and the time used for sampling}
  \item{final.values}{the final samples for each set of parameters required to
                      continue another run from the current run. Has length 5
                      if \code{algorithm = "MH"} and length 7 if \code{algorithm = "SMN"}.}
  \item{y}{the response matrix provided to the function.}
  \item{x}{the predictor matrix provided to the function.}
  \item{initial.values}{the initializing values provided to the function. If none were 
                        provided, the default initializing values generated by the
                        function.}
  \item{priors}{the priors provided to the function. If the input did not have an element
                named \code{gamma}, the default prior generated by the function
                is included as the element \code{gamma}.}
  \item{iter}{number of iterations provided to the function.}
  \item{algorithm}{the algorithm name provided to the function}
  \item{rw.scale}{the standard deviations provided to the function of the proposal distributions used for Metropolis-Hastings updates.}
  \item{debug}{logical provided to the function}
  
if \code{debug = TRUE}, the list \code{debug_out} is also included, which includes the random
variates generated during sampling.
}
\references{
  Geyer, C. J. (2011) Introduction to MCMC. In \emph{Handbook of Markov Chain Monte Carlo}. Edited 
  by S. P. Brooks, A. E. Gelman, G. L. Jones, and X. L. Meng. Chapman & Hall/CRC, Boca Raton, FL, pp. 3--48.
  
  Xiang, D. (2020) \emph{Fully Bayesian penalized regression with a generalized bridge prior}.
University of Minnesota
}
\author{
\packageAuthor{GMCB}
}
\note{
  The dimensions of \code{mu}, \code{lambda}, \code{gamma}, and \code{betas} differ depending on the 
  class of the first argument. If the first argument is a matrix, \code{mu[1,]},  
  \code{gamma[1,]}, and \code{betas[[1]]} are the initializing values passed to the function by the user,
  and \code{lambda[1,]} is a vector of 1's that serves as a placeholder because no initializing value is 
  needed. The remaining \code{iter} rows contain the samples.
  When the first argument is an object of class \code{"cov.est"}, the initializing values are taken to be 
  the final states from \code{cE}. To avoid redundancy if the user wanted to combine the output from multiple
  chains run in sequence, the initializing values are excluded, so that \code{mu[1,]}, \code{lambda[1,]}, 
  \code{gamma[1,]}, and \code{betas[[1]]} are the next set of samples, rather than the final state from the
  previous chain.
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
