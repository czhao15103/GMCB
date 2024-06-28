\name{dgbr}
\alias{dgbr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marginal density under the generalized bridge prior
}
\description{
Returns a function that can be used to plot the marginal density of a regression coefficient
under the generalized bridge prior.
}
\usage{
dgbr(e1, f1, e2, f2, k1, k2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e1}{
shape parameter of the first Gamma component of the mixture prior on the regularization parameters.
}
  \item{f1}{
rate parameter of the first Gamma component of the mixture prior on the regularization parameters.
}
  \item{e2}{
shape parameter of the second Gamma component of the mixture prior on the regularization parameters.
}
  \item{f2}{
rate parameter of the second Gamma component of the mixture prior on the regularization parameters.
}
  \item{k1}{
lower bound of the uniform prior on the penalty parameter.
}
  \item{k2}{
upper bound of the uniform prior on the penalty parameter.
}
}
\details{
The generalized bridge (GBR) prior for a regression coefficient \eqn{\beta} is as follows:
\deqn{
    \nu(\beta|\lambda, \sigma^2,\alpha) = \frac{\alpha \lambda^{1/\alpha}}{2^{1/\alpha+1}(\sigma^2)^{1/\alpha}\Gamma(1/\alpha)}e^{-(\lambda|\beta|^\alpha)/(2\sigma^2)} \\
    \lambda \sim \frac{1}{2}\text{Gamma}(e_{1},f_{1})+\frac{1}{2}\text{Gamma}(e_{2},f_{2}) \\
    \sigma^2 \sim \text{IG}(a,b) \\
    \alpha \sim \text{Unif}(k_1, k_2)
  }
This function plots the marginal density of \eqn{\beta} when \eqn{\sigma^2=1}.
}
\value{
A function that takes as its argument possible values of the regression coefficient \code{beta}.
}
\references{
    Zhao, C., Xiang, D., Jones, G. L., & Rothman, A. J. (2024). A Bayesian Generalized Bridge Regression Approach to Covariance Estimation in the Presence of Covariates. arXiv preprint arXiv:2406.00906.
}
%\author{
%Bahadir Yuzbasi, Mohammad Arashi, and Fikri Akdeniz wrote the C++ 
%source code for computing the frequentist bridge. Christina Zhao
%wrote the wrapper that adapts the frequentist bridge for 
%initializing the GMCB model.
%}
\note{
It is recommended that \code{nalpha}\eqn{\cdot}\code{nlambda} not be too large
to avoid long computation times. Even with the default values and parallelization, 
runtime can be extensive if \emph{q} is large.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}
\examples{

# lambda prior mixture of Gamma(1, 1) and Gamma(40, 0.5)
lambda.prior <- c(1, 1, 40, 0.5)

# alpha prior 
alpha.prior <- c(0.5, 4)

# produces a function
gbr1 <- dgbr(lambda.prior[1], lambda.prior[2], lambda.prior[3], lambda.prior[4],
             alpha.prior[1], alpha.prior[2])
             
x <- seq(-4, 4, length.out = 1e4)

curve(gbr1(x), type = "l", xlim = c(-1, 1), ylim = c(0, 5),
     xlab = bquote(beta), ylab = "Density")

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
