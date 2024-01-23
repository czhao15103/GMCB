\name{qda}
\alias{qda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quadratic discriminant analysis with prespecified mean and precision estimates
}
\description{
A function that performs quadratic discriminant analysis for \emph{k} classes
with user-provided class means and inverse covariance matrices.
}
\usage{
qda(test, classes, class.prior, class.means, class.precisions, parallel = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test}{
matrix of observations to be classified with dimension \emph{m} by \emph{q}.
}
  \item{classes}{
vector of length \emph{k} specifying the class labels.
}
  \item{class.prior}{
vector of length \emph{k} specifying the prior probability of each class.
Ordering is assumed to match the order of the labels provided by \code{classes}.
}
  \item{class.means}{
list of length \emph{k}, with each element a vector of length \emph{q}, specifying the class mean. 
Ordering is assumed to match the order of the labels provided by \code{classes}.
}
  \item{class.precisions}{
list of length \emph{k}, with each element a \emph{q} by \emph{q} matrix specifying the 
class precision matrix. Ordering is assumed to match the order of the labels provided by \code{classes}.
}
  \item{parallel}{
number of cores to use for parallelization. If value is less than 2,
there is no parallelization. Default value is \code{parallel = 0}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A vector of length \emph{m} with the predicted labels for each observation in \code{test}.
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

library(mvtnorm)

mats <- matgen(p = 10, rho = 0.7, type = "CS")
train <- rmvnorm(100, mean = rep(0, 10), sigma = mats$cov.mat)
test <- rmvnorm(1000, mean = rep(0, 10), sigma = mats$cov.mat)

train.class <- ifelse(train[,1] < 0, 0, 1)
test.class <- ifelse(test[,1] < 0, 0, 1)

classes <- c(0, 1)
class.prior <- c(0.5, 0.5)
class.means <- list(colMeans(train[which(train.class == 0),]),
                    colMeans(train[which(train.class == 1),]))
class.cov <- list(cov(train[which(train.class == 0),]),
                  cov(train[which(train.class == 1),]))
class.precisions <- lapply(class.cov, solve)

qda.class <- qda(test, classes, class.prior, class.means, class.precisions,
                 parallel = 0)

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
