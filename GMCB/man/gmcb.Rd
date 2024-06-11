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
gmcb(y, x = NULL, meanzero = FALSE, initial.values = NULL, priors, rw.scale, 
     iter = 1e4, algorithm = c("MH", "SMN"), debug = FALSE)
                 
\method{gmcb}{matrix}(y, x = NULL, meanzero = FALSE, initial.values = NULL, priors, rw.scale, 
     iter = 1e4, algorithm = c("MH", "SMN"), debug = FALSE)

\method{gmcb}{gbridge}(y, x = NULL, meanzero = FALSE, initial.values = NULL, priors, rw.scale, 
     iter = 1e4, algorithm = c("MH", "SMN"), debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{
response matrix. Must have at least 2 columns. 
}
  \item{x}{
optional covariate matrix, see `Details.`
}
  \item{meanzero}{
logical. Default value is \code{FALSE}. If \code{TRUE}, the argument 
\code{x} is ignored, and the model is fit assuming the mean is zero, see `Details.`
}
  \item{initial.values}{
optional list with named elements used to initialize the Markov chain, see `Details.` 
}
  \item{priors}{
list with named elements providing the necessary priors, see `Details.`
}
  \item{rw.scale}{
list with named elements giving the standard deviations of the normal proposals used 
for the Metropolis-Hastings updates, see `Details.`
}
  \item{iter}{
number of iterations. Default value is \code{1e4}.
}
  \item{debug}{
logical. Default value is \code{FALSE}. If \code{TRUE}, returns the random variates
generated during sampling as described below. 
}
  \item{algorithm}{
specifies the algorithm to run, either \code{MH} or \code{SMN}. Defaults to \code{MH}.
}
}
\details{
  Let \emph{Y} be a vector of length \emph{q}, \emph{X} a vector length \emph{p}, \emph{B} a \eqn{p \times q}
  matrix, and \eqn{\Omega} a positive-definite \eqn{q \times q} matrix. Suppose that the distribution of
  \emph{Y} is given by \eqn{N(B'X,\Omega^{-1})}, where \emph{B'} denotes the transpose of \emph{B}. 
  
  Define the modified Cholesky decomposition of \eqn{\Omega^{-1}} to be
  \deqn{
    \Omega = T'D^{-1}T
  }
  where \emph{T} is a unit lower triangular matrix and \emph{D} is a diagonal
  matrix with positive entries. By reparameterizing the likelihood of \emph{Y} using
  this modified Cholesky decomposition, the problem of estimating \eqn{(B,\Omega)}
  can be turned into estimating the regression coefficients 
  and error variances of a sequence of \emph{q} regressions.
  
  The generalized mean-covariance bridge (GMCB) prior uses the generalized bridge prior proposed
  by Xiang (2020) to model the parameters of the \emph{q} regressions. The posterior distribution
  under the GMCB prior can be sampled using either the GMCB-MH algorithm or the GMCB-SMN algorithm.
  
  The arguments \code{x} and \code{meanzero} determine the mean structure estimated by
  \code{gmcb}. While the distribution of \emph{Y} above is presented as depending on covariates,
  the function can also perform estimation for the cases \eqn{Y\sim N(\mu,\Omega^{-1})}
  and \eqn{Y\sim N(0,\Omega^{-1})}. Setting \code{meanzero = TRUE} performs estimation for
  \eqn{Y\sim N(0,\Omega^{-1})}, and the function ignores the argument \code{x}. If
  \code{meanzero = FALSE} and \code{x} is not provided, the function defaults to estimation 
  for \eqn{Y\sim N(\mu,\Omega^{-1})}. If \code{meanzero = FALSE} and \code{x} is provided,
  the function performs estimation for \eqn{Y\sim N(B'X,\Omega^{-1})}.
  
  If the argument \code{initial.values} is not provided, the function chooses initializing
  values based on the relative sizes of \emph{p}, \emph{q}, and \emph{n}. 
  Parameters with interpretations as regression coefficients are initialized by the 
  ordinary least squares (OLS) estimates when they exist and zero otherwise. Parameters 
  with interpretations as regression error variances are initialized by the OLS estimates
  when they exist and the sample variances of the columns of \code{Y} otherwise. The regularization
  parameters and penalty exponent parameters are initialized from the provided priors. 
  
  The required elements for \code{initial.values}, \code{priors}, and \code{rw.scale}
  depend on the mean structure being estimated, as well as the selected algorithm.
  When \code{meanzero = FALSE},
  \itemize{
    \item \code{priors} must have named elements \code{lambda}, \code{tau}, and \code{alpha}.
    The element \code{gamma} is optional; if excluded from the list provided to \code{gmcb},
    the function will use the method of moments to specify the prior.
    \item If \code{algorithm = MH}, 
    \itemize{
      \item \code{initial.values}, if provided, must have named elements 
      \code{b}, \code{delta}, \code{gamma}, \code{alpha.b}, and \code{alpha.d}.
      \item \code{rw.scale} must have named elements \code{b}, \code{delta}, \code{alpha.b}, 
      and \code{alpha.d}.
    }
    \item If \code{algorithm = SMN},
    \itemize{
      \item \code{initial.values}, if provided, must have named elements 
      \code{b}, \code{delta}, \code{gamma}, \code{lambda}, \code{tau}, \code{alpha.b}, 
      and \code{alpha.d}.
      \item \code{rw.scale} must have named elements \code{alpha.b} and \code{alpha.d}.
    }
  }
  When \code{meanzero = TRUE}, the required elements are as follows:
  \itemize{
    \item \code{priors} must have named elements \code{tau} and \code{alpha}.
    The element \code{gamma} is optional; if excluded from the list provided to \code{gmcb},
    the function will use the method of moments to specify the prior.
    \item If \code{algorithm = MH}, 
    \itemize{
      \item \code{initial.values}, if provided, must have named elements 
      \code{delta}, \code{gamma}, and \code{alpha.d}.
      \item \code{rw.scale} must have named elements \code{delta} and \code{alpha.d}.
    }
    \item If \code{algorithm = SMN},
    \itemize{
      \item \code{initial.values}, if provided, must have named elements 
      \code{delta}, \code{gamma}, \code{tau}, and \code{alpha.d}.
      \item \code{rw.scale} must have named element \code{alpha.d}.
    }
  }
  
  The details of the hierarchical model and the algorithms will be made available.

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
                     updates. Has length 4 if \code{algorithm = MH} and length 2
                     if \code{algorithm = SMN}.}
\item{initializing.seed}{the value of \code{.Random.seed} before generation of initializing values.}
  \item{start.seed}{the value of \code{.Random.seed} before sampling began.}
  \item{end.seed}{the value of \code{.Random.seed} after sampling ended.}
  \item{timing}{a character vector of length 2, giving the time used for 
                pre-sampling computations and the time used for sampling}
  \item{final.values}{the final samples for each set of parameters required to
                      continue another run from the current run. Has length 5
                      if \code{algorithm = MH} and length 7 if \code{algorithm = SMN}.}
  \item{y}{the response matrix provided to the function.}
  \item{x}{the predictor matrix provided to the function.}
  \item{meanzero}{logical provided to the function.}
  \item{initial.values}{the initializing values provided to the function. If none were 
                        provided, the default initializing values generated by the
                        function.}
  \item{priors}{the priors provided to the function. If the input did not have an element
                named \code{gamma}, the default prior generated by the function
                is included as the element \code{gamma}.}
  \item{iter}{number of iterations provided to the function.}
  \item{algorithm}{the algorithm name provided to the function.}
  \item{rw.scale}{the standard deviations of the proposal distributions used for 
  Metropolis-Hastings updates provided to the function.}
  \item{debug}{logical provided to the function.}
  
if \code{debug = TRUE}, the list \code{debug_out} is also included, which includes the random
variates generated during sampling. 
}
\references{
  Geyer, C. J. (2011) Introduction to MCMC. In \emph{Handbook of Markov Chain Monte Carlo}. Edited 
  by S. P. Brooks, A. E. Gelman, G. L. Jones, and X. L. Meng. Chapman & Hall/CRC, Boca Raton, FL, pp. 3--48.
  
  Pourahmadi, M. (2011) Covariance Estimation: The GLM and Regularization 
  Perspectives. \emph{Statistical Science} \bold{26}:3 369--387.
  
  Xiang, D. (2020) \emph{Fully Bayesian penalized regression with a generalized bridge prior}.
  University of Minnesota.
}
%\author{
%The code for sampling an exponentially-tilted postive stable random variable was written by Marius Hofert and Martin Maechler. The rest of the algorithm was written by Christina Zhao.
%}
\note{
  The dimensions of \code{b}, \code{delta}, \code{lambda}, \code{tau}, \code{gamma}, 
  \code{alpha.b}, and \code{alpha.d} differ depending on the class of the first argument. 
  If \code{y} is a matrix, matrix (vector) output in the component \code{mcmc} will have
  \code{iter + 1} rows (elements). The first row (element) of the matrices (vectors)
  will be the initializing values. When \code{algorithm = MH}, the first row of 
  \code{lambda} and \code{tau} will be \code{NA}, as initializing values are not required.
  
  If \code{y} has class \code{gbridge}, the initializing values are taken to be 
  the final states from \code{y}. To avoid redundancy if the user wanted to combine 
  the output a sequence of consecutive runs, the initializing values are excluded,
  so that matrix (vector) output in the component \code{mcmc} will have
  \code{iter} rows (elements).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

library(mvtnorm)

n <- 100
p <- 5
q <- 5

rho <- 0.7

x.mats <- matgen(p, rho, type = "AR")
y.mats <- matgen(q, rho, type = "AR")

x <- rmvnorm(n, mean = rep(0, p), sigma = x.mats$cov.mat)
b <- matrix(rnorm(p*q), nrow = p)
y <- x \%*\% b + rmvnorm(n, mean = rep(0, q), sigma = y.mats$cov.mat)

priors <- list(lambda = c(1, 1, 40, 0.5), tau = c(1, 1, 40, 0.5),
               alpha = c(0.5, 4))
rw.scale <- list(b = 0.1, delta = 0.1, alpha.b = 0.5, alpha.d = 0.5)

out_nodebug <- gmcb(y = y, x = x, priors = priors, rw.scale = rw.scale)
names(out_nodebug)
nrow(out_nodebug$mcmc$b)

out_debug <- gmcb(y = y, x = x, priors = priors, rw.scale = rw.scale, debug = TRUE)
names(out_debug)

## continuing a run
out2 <- gmcb(out_nodebug)
nrow(out2$mcmc$b)

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
