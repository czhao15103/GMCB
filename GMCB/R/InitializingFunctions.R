
#### GMCB initialization using frequentist bridge penalized regression ####

#### method of moments prior elicitation for gamma ####

gamma_mom <- function(gamma) {
  stopifnot(all(gamma > 0))
  
  gamma.sm <- mean(gamma)
  gamma.svar <- var(gamma)
  a <- 2 + gamma.sm^2/gamma.svar
  b <- gamma.sm + gamma.sm^3/gamma.svar
  return(c(a, b))
}

#### initialize nuisance parameters from provided priors ####

nuisance.init <- function(p, q, lambda.prior, tau.prior, alpha.prior, 
                          b.large = TRUE, delta.large = TRUE) {
  
  # error checks
  stopifnot(p > 0)
  stopifnot(q >= 2)
  stopifnot(all(lambda.prior > 0))
  stopifnot(all(tau.prior > 0))
  stopifnot(all(alpha.prior > 0))
  stopifnot(alpha.prior[1] < alpha.prior[2])
  stopifnot(is.logical(b.large))
  stopifnot(is.logical(delta.large))
  
  if ((is.null(dim(lambda.prior)) & length(lambda.prior) != 4) | !all(dim(lambda.prior) == c(4, p*q))) {
    stop("lambda.prior must either be a vector of length 4 or be a 4 by pq matrix")
  }
  
  if ((is.null(dim(tau.prior)) & length(tau.prior) != 4) | !all(dim(tau.prior) == c(4, q*(q-1)/2))) {
    stop("tau.prior must either be a vector of length 4 or be a 4 by q(q-1)/2 matrix")
  }
  
  if (is.null(dim(lambda.prior))) {
    if (b.large) { 
      # generate lambda values from first component of prior, which has small mean and variance
      
      lambda.init <- rgamma(p*q, shape = lambda.prior[1], rate = lambda.prior[2])
      
    } else {
      # generate lambda values from second component of prior, which has large mean and variance
      
      lambda.init <- rgamma(p*q, shape = lambda.prior[3], rate = lambda.prior[4])
      
    }
  } else {
    if (b.large) { 
      # generate lambda values from first component of prior, which has small mean and variance
      
      lambda.init <- rgamma(p*q, shape = lambda.prior[1,], rate = lambda.prior[2,])
      
    } else {
      # generate lambda values from second component of prior, which has large mean and variance
      
      lambda.init <- rgamma(p*q, shape = lambda.prior[3,], rate = lambda.prior[4,])
      
    }
  }
  
  if (is.null(dim(tau.prior))) {
    if (delta.large) {
      # generate tau values from first component of prior, which has small mean and variance
      
      tau.init <- rgamma(q*(q-1)/2, shape = tau.prior[1], rate = tau.prior[2])
      
    } else {
      # generate tau values from second component of prior, which has large mean and variance
      
      tau.init <- rgamma(q*(q-1)/2, shape = tau.prior[3], rate = tau.prior[4])
      
    }
  } else {
    if (delta.large) {
      # generate tau values from first component of prior, which has small mean and variance
      
      tau.init <- rgamma(q*(q-1)/2, shape = tau.prior[1,], rate = tau.prior[2,])
      
    } else {
      # generate tau values from second component of prior, which has large mean and variance
      
      tau.init <- rgamma(q*(q-1)/2, shape = tau.prior[3,], rate = tau.prior[4,])
      
    }
  }
  
  alpha.b.init <- runif(1, alpha.prior[1], alpha.prior[2])
  
  alpha.d.init <- runif(1, alpha.prior[1], alpha.prior[2])
  
  return(list(lambda.init = lambda.init, tau.init = tau.init,
              alpha.b.init = alpha.b.init, alpha.d.init = alpha.d.init))
}

nuisance.init.meanzero <- function(q, tau.prior, alpha.prior, delta.large = TRUE) {
  
  # error checks
  stopifnot(q >= 2)
  stopifnot(all(tau.prior > 0))
  stopifnot(all(alpha.prior > 0))
  stopifnot(alpha.prior[1] < alpha.prior[2])
  stopifnot(is.logical(delta.large))
  
  if ((is.null(dim(tau.prior)) & length(tau.prior) != 4) | !all(dim(tau.prior) == c(4, q*(q-1)/2))) {
    stop("tau.prior must either be a vector of length 4 or be a 4 by q(q-1)/2 matrix")
  }
  
  if (is.null(dim(tau.prior))) {
    if (delta.large) {
      # generate tau values from first component of prior, which has small mean and variance
      
      tau.init <- rgamma(q*(q-1)/2, shape = tau.prior[1], rate = tau.prior[2])
      
    } else {
      # generate tau values from second component of prior, which has large mean and variance
      
      tau.init <- rgamma(q*(q-1)/2, shape = tau.prior[3], rate = tau.prior[4])
      
    }
  } else {
    if (delta.large) {
      # generate tau values from first component of prior, which has small mean and variance
      
      tau.init <- rgamma(q*(q-1)/2, shape = tau.prior[1,], rate = tau.prior[2,])
      
    } else {
      # generate tau values from second component of prior, which has large mean and variance
      
      tau.init <- rgamma(q*(q-1)/2, shape = tau.prior[3,], rate = tau.prior[4,])
      
    }
  }
  
  alpha.d.init <- runif(1, alpha.prior[1], alpha.prior[2])
  
  return(list(tau.init = tau.init, alpha.d.init = alpha.d.init))
}

#### OLS regression initialization for just covariance parameters ####

ols.covariance <- function(y) {
  
  q <- ncol(y)
  
  # errror checks
  stopifnot(nrow(y) > q)
  
  delta.init <- numeric(q*(q-1)/2)
  gamma.init <- numeric(q)
  
  centered.y <- sweep(y, 2, colMeans(y))
  
  l <- summary(lm(centered.y[,1] ~ 1))
  gamma.init[1] <- l$sigma^2
  
  for (c in 2:q) {
    pos <- pos.vechrow(c, 1:(c-1))
    l <- summary(lm(centered.y[,c] ~ centered.y[,1:(c-1)]))
    delta.init[pos] <- l$coefficients[2:c,1]
    gamma.init[c] <- l$sigma^2
  }
  
  # gamma prior
  gamma.prior <- gamma_mom(gamma.init)
  
  return(list(delta.init = delta.init, gamma.init = gamma.init, 
              gamma.prior = gamma.prior))
}

#### OLS regression initialization ####

# using the sequence of regressions from the modified Cholesky decomposition to choose initializing
# values for the parameters
# this initialization function assumes that x has full column rank
ols.init <- function(y, x) {
  # y: n x q matrix
  # x: n x p matrix
  
  q <- ncol(y)
  lowertcount <- q*(q-1)/2 # number of strictly lower triangular elements
  
  stopifnot(q >= 2)
  stopifnot(nrow(y) == nrow(x))
  stopifnot(nrow(x) > ncol(x))
  stopifnot(nrow(y) > ncol(y))
  
  b.init.mat <- solve(crossprod(x), crossprod(x, y))  # OLS estimate
  
  # initialization
  delta.init <- numeric(lowertcount)
  gamma.init <- numeric(q)
  
  centered.y <- y - x %*% b.init.mat
  
  l <- summary(lm(centered.y[,1] ~ 1))
  gamma.init[1] <- l$sigma^2
  
  for (c in 2:q) {
    pos <- pos.vechrow(c, 1:(c-1))
    l <- summary(lm(centered.y[,c] ~ centered.y[,1:(c-1)]))
    delta.init[pos] <- l$coefficients[2:c,1]
    gamma.init[c] <- l$sigma^2
  }
  
  # gamma prior
  gamma.prior <- gamma_mom(gamma.init)
  
  return(list(b.init = as.numeric(matrixcalc::vec(b.init.mat)), 
              delta.init = delta.init, 
              gamma.init = gamma.init, gamma.prior = gamma.prior))
}

#### naive initialization ####

init_naive <- function(y, x, lambda.prior, tau.prior, alpha.prior) {
  p <- ncol(x)
  q <- ncol(y)
  
  stopifnot(nrow(y) == nrow(x))
  stopifnot(q >= 2)
  
  b.init <- rep(0, p*q)
  
  d.init <- rep(0, q*(q-1)/2)
  
  nuisance <- nuisance.init(p, q, lambda.prior, tau.prior, alpha.prior,
                            b.large = FALSE, delta.large = FALSE)
  
  gamma.init <- apply(y, 2, var)
  
  # method of moments gamma prior
  gamma.prior <- gamma_mom(gamma.init)
  
  return(list(b.init = b.init, delta.init = d.init,
              lambda.init = nuisance$lambda.init, tau.init = nuisance$tau.init,
              alpha.b.init = nuisance$alpha.b.init, alpha.d.init = nuisance$alpha.d.init,
              gamma.init = gamma.init, gamma.prior = gamma.prior))
}

#### frequentist bridge initialization ####

# # uses the C++ bridge implementation from package rbridge
#   # while the rbridge package does not allow alpha > 2, the LQA algorithm they use
#     # can be used for alpha > 2
# init_freqbridge <- function(y, x, alpha.prior, nfolds, nalpha, nlambda) {
#   
#   n <- nrow(y)
#   p <- ncol(x)
#   q <- ncol(y)
#   
#   stopifnot(nrow(x) == n)
#   stopifnot(q >= 2)
#   
#   # parameter grid values
#   alpha.grid <- seq(alpha.prior[1], alpha.prior[2], length.out = nalpha)
#   lambda.grid <- 2^((1:nlambda) - 6) # used by Park and Yoon (2011) in their simulations
#   
#   # full set of parameter values
#   pars <- expand.grid(alpha.grid, lambda.grid)
#   colnames(pars) <- c("alpha", "lambda")
#   
#   # indices for split
#   folds <- kfoldsplit(n, nfolds)
#   kfoldmse <- matrix(nrow = nfolds, ncol = nrow(pars))
#   
#   clus <- makeCluster(cores)
#   registerDoParallel(clus)
#   
#   for (k in 1:nfolds) {
#     mse <- foreach(i = 1:nrow(pars), .combine = "c", .packages = c("Matrix", "GMCB")) %dorng% {
#       
#       y.train <- y[-folds[[k]],]
#       y.test <- y[folds[[k]],]
#       x.train <- x[-folds[[k]],]
#       x.test <- x[folds[[k]],]
#       
#       fitted.b <- matrix(nrow = p, ncol = q)
#       fit.pred <- numeric(q-1)
#       
#       for (c in 2:q) {
#         y.c.train <- y.train[,c]
#         y.pred.train <- y.train[,1:(c-1),drop=FALSE] - x.train %*% fitted.b[,1:(c-1),drop=FALSE]
#         x.fit <- cbind(x.train, y.pred.train)
#         
#         y.c.test <- y.test[,c]
#         y.pred.test <- y.test[,1:(c-1),drop=FALSE] - x.test %*% fitted.b[,1:(c-1),drop=FALSE]
#         x.test <- cbind(x.test, y.pred.test)
#         
#         fit.c <- Bridge(x = x.fit, y = y.c.train, q = pars[i,1], lambda = pars[i,2], 
#                         converge = 10^10, eta = 1e-7) # default values in the rbridge package
#         fit.pred.c <- Prediction_Grid(x_test = x.test, x_train = x.fit,
#                                       y_train = y.c.train, grid_betas = as.matrix(fit.c))
#         fit.pred[c-1] <- mean((y.c.test - fit.pred.c)^2)
#       }
#       
#       return(sum(fit.pred))
#     }
#     
#     kfoldmse[k,] <- mse
#     
#   }
#   
#   stopCluster(clus)
#   
#   cvmse <- colMeans(kfoldmse)
#   cvresults <- cbind(pars, cvmse = cvmse)
#   min.par <- as.numeric(pars[which.min(cvmse),])
#   
#   # obtain the initializing values for B, delta, gamma
#   b.init <- numeric(p*q)
#   delta.init <- numeric(q*(q-1)/2)
#   gamma.init <- numeric(q)
#   gamma.init[1] <- mean((y[,1] - mean(y[,1]))^2)
#   pos <- lapply(2:q, function(c) pos.vechrow(c, 1:(c-1)))
#   
#   for (c in 2:q) {
#     x.fit <- 
#     fit.c <- Bridge(x = y[,1:(c-1),drop=FALSE], y = y[,c], q = min.par[1], lambda = min.par[2],
#                     converge = 10^10, eta = 1e-7)
#     delta.init[pos[[c-1]]] <- as.numeric(fit.c)
#     gamma.init[c] <- mean((y[,c] - y[,1:(c-1),drop=FALSE] %*% as.matrix(fit.c))^2)
#   }
#   
#   gamma.prior <- gamma_mom(gamma.init)
#   
#   out <- list(delta.init = delta.init, gamma.init = gamma.init, gamma.prior = gamma.prior,
#               alpha.b.init = min.par[1], alpha.d.init = min.par[1], 
#               lambda.init = rep(min.par[2], q), # if there are no covariates, then there are only q lambdas
#               tau.init = rep(min.par[2], q*(q-1)/2),
#               cvresults = cvresults)
#   return(out)
# }

# If there is no mean structure, the sequence of q regressions is independent
init_freqbridge_nocovariates <- function(y, alpha.prior, nfolds = 5, nalpha = 10, nlambda = 20,
                                         cores = 4, seed = NULL) {
  n <- nrow(y)
  q <- ncol(y)
  
  stopifnot(q >= 2)
  
  # parameter grid values
  alpha.grid <- seq(alpha.prior[1], alpha.prior[2], length.out = nalpha)
  lambda.grid <- 2^((1:nlambda) - 6) # used by Park and Yoon (2011) in their simulations
  
  # full set of parameter values
  pars <- expand.grid(alpha.grid, lambda.grid)
  colnames(pars) <- c("alpha", "lambda")
  
  # set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
    seedset <- .Random.seed
    assign(".Random.seed", seedset, .GlobalEnv)
  }
  
  # indices for split
  folds <- kfoldsplit(n, nfolds)
  kfoldmse <- matrix(nrow = nfolds, ncol = nrow(pars))
  
  clus <- makeCluster(cores)
  registerDoParallel(clus)
  
  for (k in 1:nfolds) {
    mse <- foreach(i = 1:nrow(pars), .combine = "c", .packages = c("Matrix", "GMCB")) %dorng% {
      
      train <- y[-folds[[k]],]
      test <- y[folds[[k]],]
      
      fit.pred <- numeric(q-1)
      
      for (c in 2:q) {
        y.c.train <- train[,c]
        x.c.train <- train[,1:(c-1),drop=FALSE]
        y.c.test <- test[,c]
        x.c.test <- test[,1:(c-1),drop=FALSE]
        
        fit.c <- Bridge(x = x.c.train, y = y.c.train, q = pars[i,1], lambda = pars[i,2], 
                        converge = 10^10, eta = 1e-7) # default values in the rbridge package
        fit.pred.c <- Prediction_Grid(x_test = x.c.test, x_train = x.c.train,
                                      y_train = y.c.train, grid_betas = as.matrix(fit.c))
        fit.pred[c-1] <- mean((y.c.test - fit.pred.c)^2)
      }
      
      return(sum(fit.pred))
    }
    
    kfoldmse[k,] <- mse
    
  }
  
  stopCluster(clus)
  
  cvmse <- colMeans(kfoldmse)
  cvresults <- cbind(pars, cvmse = cvmse)
  min.par <- as.numeric(pars[which.min(cvmse),])
  
  # obtain the initializing values for delta, gamma
  delta.init <- numeric(q*(q-1)/2)
  gamma.init <- numeric(q)
  gamma.init[1] <- mean((y[,1] - mean(y[,1]))^2)
  pos <- lapply(2:q, function(c) pos.vechrow(c, 1:(c-1)))
  
  for (c in 2:q) {
    fit.c <- Bridge(x = y[,1:(c-1),drop=FALSE], y = y[,c], q = min.par[1], lambda = min.par[2],
                    converge = 10^10, eta = 1e-7)
    delta.init[pos[[c-1]]] <- as.numeric(fit.c)
    gamma.init[c] <- mean((y[,c] - y[,1:(c-1),drop=FALSE] %*% as.matrix(fit.c))^2)
  }
  
  gamma.prior <- gamma_mom(gamma.init)
  
  out <- list(delta.init = delta.init, gamma.init = gamma.init, gamma.prior = gamma.prior,
              alpha.b.init = min.par[1], alpha.d.init = min.par[1], 
              lambda.init = rep(min.par[2], q), # if there are no covariates, then there are only q lambdas
              tau.init = rep(min.par[2], q*(q-1)/2),
              cvresults = cvresults)
  return(out)
}
