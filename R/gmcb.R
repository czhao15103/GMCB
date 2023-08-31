
#### GMCB wrapper function ####

gmcb <- function(y, x = NULL, 
                 initial.values = NULL,
                 priors,
                 rw.scale,
                 iter = 1e4, 
                 algorithm = c("MH", "SMN"),
                 debug = FALSE) {
  UseMethod("gmcb")
}

#### matrix method ####

gmcb.matrix <- function(y, x = NULL, 
                        initial.values = NULL,
                        priors,
                        rw.scale,
                        iter = 1e4, 
                        algorithm = c("MH", "SMN"),
                        debug = FALSE) {
  
  #### error checking for data and algorithm choice ####
  
  if (!(algorithm %in% c("MH", "SMN"))) {
    stop("algorithm must be either 'MH' or 'SMN'")
  }
  
  if (!("matrix" %in% class(y)))
    stop("y must be a matrix")
  
  n <- nrow(y)
  q <- ncol(y)
  
  if (q < 2) {
    stop("y must have at least 2 columns")
  }
  
  if (is.null(x)) { # if no covariates supplied, will perform intercept-only regression
    x <- matrix(1, nrow = n, ncol = 1)
  } else if (is.null(dim(x))) { # if x is provided as a numeric vector
    x <- matrix(x, ncol = 1)
  } 
  
  p <- ncol(x)
  pq <- p*q
  lowertcount <- q*(q-1)/2 # number of strictly lower triangular elements in a q x q matrix
  iter.1 <- iter + 1 # will include initial values as "step 0"
  
  if (nrow(x) != n) stop("x and y must have the same number of rows")
  
  ##### error checks for priors #####
  required.priors.names <- c("lambda", "tau", "alpha") # gamma is optional
  names.in.prior.input <- required.priors.names %in% names(priors)
  if (!all(names.in.prior.input)) {
    stop(paste("priors must include elements named", 
               paste(required.priors.names[!names.in.prior.input], collapse = ", ")))
  }
  
  # check that elements of priors are the correct class
  priors.class <- lapply(priors, class)
  list.in.priors.class <- lapply(priors.class, grepl, pattern = "list", fixed = TRUE)
  if (any(unlist(list.in.priors.class)))
    stop("elements of priors cannot be lists. They should either be numeric vectors or matrices")
  
  # alpha prior
  if (!is.null(dim(priors$alpha))) 
    stop("priors$alpha should be a numeric vector, not a matrix")
  if (length(priors$alpha) != 2) 
    stop("priors$alpha should be vector of length 2")
  if (any(priors$alpha < 0)) 
    stop("Values of priors$alpha must be nonnegative")
  if (priors$alpha[1] > priors$alpha[2]) {
    warning("First element of priors$alpha should be the lower bound. Swapping order of values.")
    priors$alpha <- rev(priors$alpha) # if first element greater than second element, switch the order
  }
  
  # lambda prior
  if ((is.null(dim(priors$lambda)) & length(priors$lambda) != 4) | !all(dim(priors$lambda) == c(4, pq))) {
    stop("priors$lambda must either be a vector of length 4 or be a 4 by ncol(x) * ncol(y) matrix")
  }
  if (!all(priors$lambda > 0)) 
    stop("All values in priors$lambda must be greater than zero")
  
  # tau prior
  if ((is.null(dim(priors$tau)) & length(priors$tau) != 4) | !all(dim(priors$tau) == c(4, lowertcount))) {
    stop("priors$tau must either be a vector of length 4 or be a 4 by ncol(y)(ncol(y) - 1) matrix")
  }
  if (!all(priors$tau > 0)) 
    stop("All values in priors$tau must be greater than zero")
  
  #### behavior if no initializing values provided ####
  
  if (is.null(initial.values)) {
    if (p < n & q < n) {
      # initialize with OLS estimates for B and delta
      
      init <- ols.init(y, x)
      
      nuisance <- nuisance.init(p, q, priors$lambda, priors$tau, priors$alpha,
                                b.large = TRUE, delta.large = TRUE)
      
      initial.values <- list(b = init$b.init, delta = init$delta.init,
                             gamma = init$gamma.init,
                             alpha.b = nuisance$alpha.b.init,
                             alpha.d = nuisance$alpha.d.init)
      
      if (algorithm == "SMN") {
        initial.values <- c(initial.values, list(lambda = nuisance$lambda.init,
                                                 tau = nuisance$tau.init))
      }
      
    } else if (p < n) { # must be the case that q > n
      # initialize with OLS estimate for B and naive initialization for covariance
      
      b.init <- solve(crossprod(x), crossprod(x, y))
      
      delta.init <- rep(0, lowertcount)
      
      nuisance <- nuisance.init(p, q, priors$lambda, priors$tau, priors$alpha,
                                b.large = TRUE, delta.large = FALSE)
      
      initial.values <- list(b = as.numeric(matrixcalc::vec(b.init)), delta = delta.init,
                             gamma = apply(y, 2, var),
                             alpha.b = nuisance$alpha.b.init,
                             alpha.d = nuisance$alpha.d.init)
      
      if (algorithm == "SMN") {
        initial.values <- c(initial.values, list(lambda = nuisance$lambda.init,
                                                 tau = nuisance$tau.init))
      }
      
    } else if (q < n) { # must be the case that p > n
      # naive initialization for B and OLS for delta
      
      b.init <- rep(0, pq)
      
      covariance.init <- ols.covariance(y)
      
      nuisance <- nuisance.init(p, q, priors$lambda, priors$tau, priors$alpha,
                                b.large = FALSE, delta.large = TRUE)
      
      initial.values <- list(b = b.init, delta = covariance.init$delta.init,
                             gamma = covariance.init$gamma.init,
                             alpha.b = nuisance$alpha.b.init,
                             alpha.d = nuisance$alpha.d.init)
      
      if (algorithm == "SMN") {
        initial.values <- c(initial.values, list(lambda = nuisance$lambda.init,
                                                 tau = nuisance$tau.init))
      }
      
    } else {
      # naive initialization for both B and delta
      
      init <- init_naive(y, x, priors$lambda, priors$tau, priors$alpha)
      
      initial.values <- list(b = init$b.init, delta = init$delta.init,
                             gamma = init$gamma.init,
                             alpha.b = init$alpha.b.init,
                             alpha.d = init$alpha.d.init)
      
      if (algorithm == "SMN") {
        initial.values <- c(initial.values, list(lambda = init$lambda.init,
                                                 tau = init$tau.init))
      }
    }
    
  }
  
  ##### error checks for initializing values and random walk scale values #####
  initial.values.names <- c("b", "delta", "gamma", "lambda", "tau", "alpha.b", "alpha.d")
  names.in.initial.values.input <- initial.values.names %in% names(initial.values)
  
  rw.scale.names <- c("b", "delta", "alpha.b", "alpha.d")
  names.in.rw.scale.input <- rw.scale.names %in% names(rw.scale)
  
  if (algorithm == "MH") {
    # GMCB-MH does not require initializing values for lambda and tau
    # requires scale for rw updates for B, delta, alpha.b, and alpha.d
    
    if (!all(names.in.initial.values.input[-c(4:5)])) { 
      
      stop(paste("When algorithm = MH, initial.values must include elements named", 
                 paste(initial.values.names[-c(4:5)][!names.in.initial.values.input[-c(4:5)]], 
                       collapse = ", ")))
    }
    
    if (!all(names.in.rw.scale.input)) { 
      
      stop(paste("When algorithm = MH, rw.scale must include elements named", 
                 paste(rw.scale.names[!names.in.rw.scale.input], collapse = ", ")))
    }
    
    if (length(rw.scale$b) == 1) {
      rw.scale$b <- rep(rw.scale$b, pq)
    }
    
    if (length(rw.scale$b) != pq) 
      stop("rw.scale$b must be a vector of length ncol(x) * ncol(y)")
    
    if (!all(rw.scale$b > 0)) {
      stop("All values of rw.scale$b must be positive")
    }
    
    if (length(rw.scale$delta) == 1) {
      rw.scale$delta <- rep(rw.scale$delta, lowertcount)
    }
    
    if (length(rw.scale$delta) != lowertcount) 
      stop("rw.scale$delta must be a vector of length ncol(y) * (ncol(y) - 1)/2")
    
    if (!all(rw.scale$delta > 0)) {
      stop("All values of rw.scale$delta must be positive")
    }
    
  } else {
    # GMCB-SMN requires initializing values for all parameters
    # requires scale for rw updates for alpha.b, and alpha.d
    
    if (!all(names.in.initial.values.input)) {
      stop(paste("When algorithm = SMN, initial.values must include elements named", 
                 paste(initial.values.names[!names.in.initial.values.input], collapse = ", ")))
    }
    
    if (length(initial.values$lambda) == 1) 
      initial.values$lambda <- rep(initial.values$lambda, pq)
    
    if (length(initial.values$lambda) != pq) 
      stop("initial.values$lambda must be a vector of length ncol(x) * ncol(y)")
    
    if (!all(initial.values$lambda > 0)) 
      stop("all values in initial.values$lambda must be greater than zero")
    
    if (length(initial.values$tau) == 1) 
      initial.values$tau <- rep(initial.values$tau, lowertcount)
    
    if (length(initial.values$tau) != lowertcount) 
      stop("initial.values$tau must be a vector of length ncol(y) * (ncol(y) - 1)/2")
    
    if (!all(initial.values$tau > 0)) 
      stop("all values in initial.values$tau must be greater than zero")
    
    if (!all(names.in.rw.scale.input[3:4])) { 
      
      stop(paste("When algorithm = SMN, initial.values must include elements named", 
                 paste(rw.scale.names[3:4][!names.in.rw.scale.input[3:4]], collapse = ", ")))
    }
  }
  
  ## common initializing values ##
  if (length(initial.values$b) != pq | !is.null(dim(initial.values$b))) 
    stop("initial.values$b must a vector of length ncol(x) * ncol(y)")
  
  if (length(initial.values$delta) != lowertcount | !is.null(dim(initial.values$delta))) 
    stop("initial.values$delta must be a vector of length ncol(y) * (ncol(y) - 1)/2")
  
  if (length(initial.values$gamma) != q | !is.null(dim(initial.values$gamma))) 
    stop("initial.values$gamma must be a vector of length ncol(y)")
  if (any(initial.values$gamma <= 0)) 
    stop("initial.values$gamma values must be greater than 0")
  
  if (length(initial.values$alpha.b) != 1) 
    stop("initial.values$alpha.b must be a scalar")
  if (initial.values$alpha.b <= 0) 
    stop("Must have initial.values$alpha.b > 0")
  
  if (length(initial.values$alpha.d) != 1) 
    stop("initial.values$alpha.d must be a scalar")
  if (initial.values$alpha.d <= 0) 
    stop("Must have initial.values$alpha.d > 0")
  
  if (initial.values$alpha.b < priors$alpha[1] | initial.values$alpha.b > priors$alpha[2])
    stop("invalid value for initial.values$alpha.b; must be between first and second elements of priors$alpha")
  
  if (initial.values$alpha.d < priors$alpha[1] | initial.values$alpha.d > priors$alpha[2])
    stop("invalid value for initial.values$alpha.d; must be between first and second elements of priors$alpha")
  
  ## gamma prior - does not have to be provided
  if (!("gamma" %in% names(priors))) {
    priors$gamma <- gamma.mom(initial.values$gamma)
  }
  
  if (!is.null(dim(priors$gamma)))
    stop("priors$gamma should be a numeric vector, not a matrix")
  if (length(priors$gamma) != 2) 
    stop("priors$gamma should be vector of length 2")
  if (any(priors$gamma <= 0)) 
    stop("All values in priors$gamma must be greater than zero")
  
  ## common rw scale ##
  if (length(rw.scale$alpha.b) != 1) 
    stop("rw.scale$alpha.b must be a scalar")
  if (rw.scale$alpha.b <= 0) 
    stop("Must have rw.scale$alpha.b > 0")
  
  if (length(rw.scale$alpha.d) != 1) 
    stop("rw.scale$alpha.d must be a scalar")
  if (rw.scale$alpha.d <= 0) 
    stop("Must have rw.scale$alpha.d > 0")
  
  #### reorganizing arguments as needed ####
  if (is.null(dim(priors$lambda))) {
    priors$lambda <- matrix(rep(priors$lambda, pq), nrow = 4, ncol = pq)
  }
  
  if (is.null(dim(priors$tau))) {
    priors$tau <- matrix(rep(priors$tau, lowertcount), nrow = 4, ncol = lowertcount)
  }
  
  #### preparation ####
  pos <- lapply(2:q, function(c) pos.vechrow(c, 1:(c-1))) # indexing for each delta_c, list of length q-1
  
  # create .Random.seed if it does not exist
  if (!exists(".Random.seed")) {
    rnorm(1)
  }
  starting.seed <- .Random.seed
  
  if (algorithm == "MH") {
    
    if (debug) {
      out <- gmcb_mh_debug(y = y, x = x, 
                           b_init = initial.values$b, 
                           d_init = initial.values$delta, 
                           gamma_init = initial.values$gamma, 
                           alpha_b_init = initial.values$alpha.b, 
                           alpha_d_init = initial.values$alpha.d, 
                           gamma_prior = priors$gamma, 
                           alpha_prior = priors$alpha, 
                           lambda_prior = priors$lambda,
                           tau_prior = priors$tau, 
                           iter = iter, 
                           b_scale = rw.scale$b, 
                           d_scale = rw.scale$delta,
                           alpha_b_scale = rw.scale$alpha.b,
                           alpha_d_scale = rw.scale$alpha.d, 
                           pos = pos)
    } else {
      out <- gmcb_mh_nodebug(y = y, x = x, 
                             b_init = initial.values$b, 
                             d_init = initial.values$delta, 
                             gamma_init = initial.values$gamma, 
                             alpha_b_init = initial.values$alpha.b, 
                             alpha_d_init = initial.values$alpha.d, 
                             gamma_prior = priors$gamma, 
                             alpha_prior = priors$alpha, 
                             lambda_prior = priors$lambda,
                             tau_prior = priors$tau, 
                             iter = iter, 
                             b_scale = rw.scale$b, 
                             d_scale = rw.scale$delta,
                             alpha_b_scale = rw.scale$alpha.b,
                             alpha_d_scale = rw.scale$alpha.d, 
                             pos = pos)
    }
    
  } else {
    
    if (debug) {
      out <- gmcb_smn_debug(y = y, x = x, 
                            b_init = initial.values$b, 
                            d_init = initial.values$delta, 
                            gamma_init = initial.values$gamma, 
                            alpha_b_init = initial.values$alpha.b, 
                            alpha_d_init = initial.values$alpha.d, 
                            lambda_init = initial.values$lambda,
                            tau_init = initial.values$tau,
                            gamma_prior = priors$gamma, 
                            alpha_prior = priors$alpha, 
                            lambda_prior = priors$lambda,
                            tau_prior = priors$tau, 
                            iter = iter, 
                            alpha_b_scale = rw.scale$alpha.b,
                            alpha_d_scale = rw.scale$alpha.d, 
                            pos = pos)
    } else {
      out <- gmcb_smn_nodebug(y = y, x = x, 
                              b_init = initial.values$b, 
                              d_init = initial.values$delta, 
                              gamma_init = initial.values$gamma, 
                              alpha_b_init = initial.values$alpha.b, 
                              alpha_d_init = initial.values$alpha.d, 
                              lambda_init = initial.values$lambda,
                              tau_init = initial.values$tau,
                              gamma_prior = priors$gamma, 
                              alpha_prior = priors$alpha, 
                              lambda_prior = priors$lambda,
                              tau_prior = priors$tau, 
                              iter = iter, 
                              alpha_b_scale = rw.scale$alpha.b,
                              alpha_d_scale = rw.scale$alpha.d, 
                              pos = pos)
    }
    
  }
  
  ending.seed <- .Random.seed
  
  ## format timing output - returned from C++ as seconds
  seconds.per.step <- diff(out$timing) # will always have length 3
  units.per.step <- numeric(2)
  for (t in 1:2) {
    if (seconds.per.step[t] < 60) {
      units.per.step[t] <- " secs"
    } else if (seconds.per.step[t] < 3600) {
      units.per.step[t] <- " mins"
      seconds.per.step[t] <- seconds.per.step[t]/60 # convert to minutes
    } else if (seconds.per.step[t] < 86400) {
      units.per.step[t] <- " hrs"
      seconds.per.step[t] <- seconds.per.step[t]/3600 # convert to hours
    } else {
      units.per.step[t] <- " days"
      seconds.per.step[t] <- seconds.per.step[t]/3600 # convert to days
    }
  }
  timing <- paste0(round(seconds.per.step, digits = 4), units.per.step)
  names(timing) <- c("Pre-sampling computations", "Sampling")
  
  ## format starting values for continuing runs
  final.values <- list(b = out$mcmc$b[iter.1,], delta = out$mcmc$delta[iter.1,],
                       gamma = out$mcmc$gamma[iter.1,], alpha.b = out$mcmc$alpha_b[iter.1],
                       alpha.d = out$mcmc$alpha_d[iter.1])
  
  if (algorithm == "SMN") {
    final.values <- c(final.values, list(lambda = out$mcmc$lambda[iter.1,],
                      tau = out$mcmc$tau[iter.1,]))
    acceptances = out$acceptances
  } else {
    acceptances <- list(b_accept_rate = colMeans(out$acceptances$b_accept_rate),
                        delta_accept_rate = colMeans(out$acceptances$delta_accept_rate),
                        alpha_b_accept_rate = mean(out$acceptances$alpha_b_accept_rate),
                        alpha_d_accept_rate = mean(out$acceptances$alpha_d_accept_rate))
  }
  
  final.out <- list(mcmc = out$mcmc, acceptances = acceptances,
                    start.seed = starting.seed, end.seed = ending.seed,
                    timing = timing, final.values = final.values,
                    y = y, x = x, initial.values = initial.values,
                    priors = priors, iter = iter, algorithm = algorithm, 
                    rw.scale = rw.scale, debug = debug)
  
  if (debug) {
    final.out <- c(final.out, list(debug_out = out$debug_out))
  }
  
  attr(final.out, "class") <- "gbridge"
  return(final.out)
}

#### gbridge method ####

# pick up run where the previous run left off

gmcb.gbridge <- function(y, x = NULL, 
                         initial.values = NULL,
                         priors,
                         rw.scale,
                         iter = 1e4, 
                         algorithm = c("MH", "SMN"),
                         debug = FALSE) {
  
  if (missing(priors)) priors <- y$priors
  if (missing(rw.scale)) rw.scale <- y$rw.scale
  if (missing(iter)) iter <- y$iter
  if (missing(algorithm)) algorithm <- y$algorithm
  if (missing(debug)) debug <- y$debug
  
  response <- y$y
  x <- y$x
  
  initial.values <- y$final.values # final values from previous run are now initializing values
  
  assign(".Random.seed", y$end.seed, .GlobalEnv)
  
  out <- gmcb.matrix(response, x = x, 
                     initial.values = initial.values,
                     priors = priors,
                     rw.scale = rw.scale,
                     iter = iter, 
                     algorithm = algorithm,
                     debug = debug)
  
  # since initialized using final values from previous run, remove initial values from matrix output 
  out$mcmc$b <- out$mcmc$b[-1,,drop = FALSE]
  out$mcmc$delta <- out$mcmc$delta[-1,,drop = FALSE]
  out$mcmc$gamma <- out$mcmc$gamma[-1,,drop = FALSE]
  out$mcmc$lambda <- out$mcmc$lambda[-1,,drop = FALSE]
  out$mcmc$tau <- out$mcmc$tau[-1,,drop = FALSE]
  out$mcmc$alpha_b <- out$mcmc$alpha_b[-1]
  out$mcmc$alpha_d <- out$mcmc$alpha_d[-1]
  return(out)
}
  