
#### GMCB wrapper function ####

gmcb <- function(y, x = NULL, meanzero = FALSE,
                 initial.values = NULL,
                 priors,
                 rw.scale,
                 iter = 1e4, 
                 algorithm = c("MH", "SMN"),
                 debug = FALSE) {
  UseMethod("gmcb")
}

#### matrix method ####

gmcb.matrix <- function(y, x = NULL, meanzero = FALSE,
                        initial.values = NULL,
                        priors,
                        rw.scale,
                        iter = 1e4, 
                        algorithm = c("MH", "SMN"),
                        debug = FALSE) {
  
  #### error checking for response and algorithm choice ####
  
  if (!(algorithm %in% c("MH", "SMN"))) {
    stop("algorithm must be either 'MH' or 'SMN'")
  }
  
  if (!("matrix" %in% class(y)))
    stop("y must be a matrix")
  
  if (length(meanzero) > 1) {
    meanzero <- meanzero[1]
  }
  
  if (length(iter) > 1) {
    iter <- iter[1]
  } 
  
  if (length(debug) > 1) {
    debug <- debug[1]
  }
  
  stopifnot(is.logical(meanzero), is.logical(debug), is.numeric(iter))
  
  n <- nrow(y)
  q <- ncol(y)
  
  if (q < 2) {
    stop("y must have at least 2 columns")
  }
  
  lowertcount <- q*(q-1)/2 # number of strictly lower triangular elements in a q x q matrix
  
  #### check that priors and rw.scale are numeric matrices/vectors ####
  
  # priors
  priors.wrong.class <- sapply(priors, is.list)
  if (any(priors.wrong.class))
    stop("elements of priors should either be vectors or matrices")
  
  priors.not.numeric <- sapply(priors, function(x) !is.numeric(x))
  if (any(priors.not.numeric))
    stop("elements of priors should be numeric")
  
  # rw.scale
  rw.scale.wrong.class <- sapply(rw.scale, is.list)
  if (any(rw.scale.wrong.class))
    stop("elements of rw.scale should be vectors")
  
  rw.scale.not.numeric <- sapply(rw.scale, function(x) !is.numeric(x))
  if (any(rw.scale.not.numeric))
    stop("elements of rw.scale should be numeric")
  
  #### checking quantities unique to nonzero mean that apply for both algorithms ####
  
  if (!meanzero) { # if meanzero = FALSE
    
    ##### checking x #####
    if (is.null(x)) { # if no covariates supplied, will perform intercept-only regression
      x <- matrix(1, nrow = n, ncol = 1)
    } else if (is.null(dim(x))) { # if x is provided as a numeric vector
      x <- matrix(x, ncol = 1)
    } 
    
    p <- ncol(x)
    pq <- p*q
    
    if (nrow(x) != n) stop("x and y must have the same number of rows")
    
    ##### lambda prior #####
    
    # checking this prior is only necessary if there mean is nonzero
    
    # check present
    if (!("lambda" %in% names(priors))) {
      stop("When meanzero = FALSE, priors must include element named lambda")
    }
    
    # check input format 
    if ((is.null(dim(priors$lambda)) & length(priors$lambda) != 4) | !all(dim(priors$lambda) == c(4, pq))) {
      stop("priors$lambda must either be a vector of length 4 or be a 4 by ncol(x) * ncol(y) matrix")
    }
    
    # check validity
    if (!all(priors$lambda > 0)) 
      stop("All values in priors$lambda must be greater than zero")
    
    # restructure if needed
    if (is.null(dim(priors$lambda))) {
      priors$lambda <- matrix(rep(priors$lambda, pq), nrow = 4, ncol = pq)
    }
    
    ##### scaling for b #####
    
    # only necessary if algorithm = "MH"
    
    # check presence
    if (algorithm == "MH") {
      
      # check presence 
      if (!("b" %in% names(rw.scale))) {
        stop("When meanzero = FALSE and algorithm = MH, rw.scale must include element named b")
      }
      
      # check formatting
      if (length(rw.scale$b) == 1) {
        rw.scale$b <- rep(rw.scale$b, pq)
      }
      
      if (length(rw.scale$b) != pq) 
        stop("rw.scale$b must be a vector of length ncol(x) * ncol(y)")
      
      # check validity
      if (!all(rw.scale$b > 0)) {
        stop("All values of rw.scale$b must be positive")
      }
    }
    
    ##### scaling for alpha.b #####
    
    # check presence
    if (!("alpha.b" %in% names(rw.scale))) {
      stop("When meanzero = FALSE, rw.scale must include element named alpha.b")
    }
    
    # check format
    if (length(rw.scale$alpha.b) < 1) {
      stop("rw.scale$alpha.b must be a scalar")
    } else if (length(rw.scale$alpha.b) > 1) {
      rw.scale$alpha.b <- rw.scale$alpha.b[1] # ignore all but first element if more than value provided
    }
    
    # check validity
    if (rw.scale$alpha.b <= 0) 
      stop("Must have rw.scale$alpha.b > 0")
  } 
  
  #### checking remaining priors ####
  
  ##### tau prior #####
  
  # check present
  if (!("tau" %in% names(priors))) {
    stop("priors must include element named tau")
  }
  
  # check input format
  if ((is.null(dim(priors$tau)) & length(priors$tau) != 4) | !all(dim(priors$tau) == c(4, lowertcount))) {
    stop("priors$tau must either be a vector of length 4 or be a 4 by ncol(y)(ncol(y) - 1) matrix")
  }
  
  # check validity
  if (!all(priors$tau > 0)) 
    stop("All values in priors$tau must be greater than zero")
  
  # restructure if needed 
  if (is.null(dim(priors$tau))) {
    priors$tau <- matrix(rep(priors$tau, lowertcount), nrow = 4, ncol = lowertcount)
  }
  
  ###### alpha prior ######
  
  # check present
  if (!("alpha" %in% names(priors))) {
    stop("priors must include element named alpha")
  }
  
  # check structure
  if (!is.null(dim(priors$alpha))) 
    stop("priors$alpha should be a numeric vector, not a matrix")
  if (length(priors$alpha) != 2) 
    stop("priors$alpha should be vector of length 2")
  
  # check validity
  if (any(priors$alpha < 0)) 
    stop("Values of priors$alpha must be nonnegative")
  if (priors$alpha[1] > priors$alpha[2]) {
    warning("First element of priors$alpha should be the lower bound. Swapping order of values.",
            immediate. = TRUE)
    priors$alpha <- rev(priors$alpha) # if first element greater than second element, switch the order
  }
  
  if (algorithm == "SMN" & priors$alpha[2] > 2) {
    stop("When algorithm = SMN, cannot have priors$alpha[2] > 2")
  }
  
  #### remaining rw.scale checks ####
  
  ##### alpha.d #####
  
  # check presence
  if (!("alpha.d" %in% names(rw.scale))) {
    stop("rw.scale must include element named alpha.d")
  }
  
  # check format
  if (length(rw.scale$alpha.d) < 1) {
    stop("rw.scale$alpha.d must be a scalar")
  } else if (length(rw.scale$alpha.d) > 1) {
    rw.scale$alpha.d <- rw.scale$alpha.d[1] # ignore all but first element if more than value provided
  }

  # check validity 
  if (rw.scale$alpha.d <= 0) 
    stop("Must have rw.scale$alpha.d > 0")
  
  ##### delta #####
  
  # only necessary if algorithm = MH
  
  if (algorithm == "MH") {
    
    # check presence
    if (!("delta" %in% names(rw.scale))) {
      stop("When algorithm = MH, rw.scale must include element named delta")
    }
    
    # check format
    if (length(rw.scale$delta) == 1) {
      rw.scale$delta <- rep(rw.scale$delta, lowertcount)
    }
    
    if (length(rw.scale$delta) != lowertcount) 
      stop("rw.scale$delta must be a vector of length ncol(y) * (ncol(y) - 1)/2")
    
    # check validity
    if (!all(rw.scale$delta > 0)) {
      stop("All values of rw.scale$delta must be positive")
    }
  }
  
  #### initialize if needed ####
  
  if (is.null(initial.values) & !meanzero) { # initialization for nonzero mean structure
    if (p < n & q < n) {
      # initialize with OLS estimates for B and delta
      
      init <- ols.init(y, x)
      
      nuisance <- nuisance.init(p, q, priors$lambda, priors$tau, priors$alpha,
                                b.large = TRUE, delta.large = TRUE)
      
      initial.values <- list(b = init$b.init, 
                             delta = init$delta.init,
                             gamma = init$gamma.init,
                             alpha.b = nuisance$alpha.b.init,
                             alpha.d = nuisance$alpha.d.init)
      
      if (algorithm == "SMN") {
        initial.values <- c(initial.values, list(lambda = nuisance$lambda.init,
                                                 tau = nuisance$tau.init))
      }
      
    } else if (p < n) { # must be the case that q >= n
      # initialize with OLS estimate for B and naive initialization for covariance
      
      b.init <- solve(crossprod(x), crossprod(x, y))
      
      delta.init <- rep(0, lowertcount)
      
      nuisance <- nuisance.init(p, q, priors$lambda, priors$tau, priors$alpha,
                                b.large = TRUE, delta.large = FALSE)
      
      initial.values <- list(b = as.numeric(matrixcalc::vec(b.init)), 
                             delta = delta.init,
                             gamma = apply(y, 2, var),
                             alpha.b = nuisance$alpha.b.init,
                             alpha.d = nuisance$alpha.d.init)
      
      if (algorithm == "SMN") {
        initial.values <- c(initial.values, list(lambda = nuisance$lambda.init,
                                                 tau = nuisance$tau.init))
      }
      
    } else if (q < n) { # must be the case that p >= n
      # naive initialization for B and OLS for delta
      
      b.init <- rep(0, pq)
      
      covariance.init <- ols.covariance(y)
      
      nuisance <- nuisance.init(p, q, priors$lambda, priors$tau, priors$alpha,
                                b.large = FALSE, delta.large = TRUE)
      
      initial.values <- list(b = b.init, 
                             delta = covariance.init$delta.init,
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
      
      initial.values <- list(b = init$b.init, 
                             delta = init$delta.init,
                             gamma = init$gamma.init,
                             alpha.b = init$alpha.b.init,
                             alpha.d = init$alpha.d.init)
      
      if (algorithm == "SMN") {
        initial.values <- c(initial.values, list(lambda = init$lambda.init,
                                                 tau = init$tau.init))
      }
    }
  } else if (is.null(initial.values) & meanzero) { # initialization when mean = 0
    if (q < n) {
      covariance.init <- ols.covariance(y)
      
      nuisance <- nuisance.init.meanzero(q, priors$tau, priors$alpha, delta.large = TRUE)
      
      initial.values <- list(delta = covariance.init$delta.init,
                             gamma = covariance.init$gamma.init,
                             alpha.d = nuisance$alpha.d.init)
      
      if (algorithm == "SMN") {
        initial.values <- c(initial.values, list(tau = nuisance$tau.init))
      }
    } else {
      delta.init <- rep(0, lowertcount)
      
      nuisance <- nuisance.init.meanzero(q, priors$tau, priors$alpha, delta.large = TRUE)
      
      initial.values <- list(delta = delta.init,
                             gamma = apply(y, 2, var),
                             alpha.d = nuisance$alpha.d.init)
      
      if (algorithm == "SMN") {
        initial.values <- c(initial.values, list(tau = nuisance$tau.init))
      }
    }
  }
  
  #### checking initial.values #####
  
  # check that they are numeric vectors
  initial.values.wrong.class <- sapply(initial.values, is.list)
  if (any(initial.values.wrong.class))
    stop("elements of initial.values should be vectors")
  
  initial.values.not.numeric <- sapply(initial.values, function(x) !is.numeric(x))
  if (any(initial.values.not.numeric))
    stop("elements of initial.values should be numeric")
  
  ##### nonzero mean #####
  
  # check b, lambda, alpha.b initial values
  
  if (!meanzero) {
    ###### b ######
    
    # check presence
    if (!("b" %in% names(initial.values))) {
      stop("When meanzero = FALSE, initial.values must have element named b")
    }
    
    # check format
    if (length(initial.values$b) == 1) {
      initial.values$b <- rep(initial.values$b, pq)
    }
    
    if (length(initial.values$b) != pq | !is.null(dim(initial.values$b))) 
      stop("initial.values$b must a vector of length ncol(x) * ncol(y)")
    
    ###### lambda ######
    
    # only required if algorithm = SMN
    
    if (algorithm == "SMN") {
      
      # check presence
      if (!("lambda" %in% names(initial.values))) {
        stop("When meanzero = FALSE and algorithm = SMN, initial.values must have element named lambda")
      }
      
      # check format
      if (length(initial.values$lambda) == 1) 
        initial.values$lambda <- rep(initial.values$lambda, pq)
      
      if (length(initial.values$lambda) != pq) 
        stop("initial.values$lambda must be a vector of length ncol(x) * ncol(y)")
      
      # check validity
      if (!all(initial.values$lambda > 0)) 
        stop("all values in initial.values$lambda must be greater than zero")
    }
    
    ###### alpha.b ######
    
    # check presence
    if (!("alpha.b" %in% names(initial.values))) {
      stop("When meanzero = FALSE, initial.values must have element named alpha.b")
    }
    
    # check format
    if (length(initial.values$alpha.b) < 1) {
      stop("initial.values$alpha.b must be a scalar")
    } else if (length(initial.values$alpha.b) > 1) {
      initial.values$alpha.b <- initial.values$alpha.b[1] # ignore all but first element if more than value provided
    }
      
    # check validity
    if (initial.values$alpha.b <= 0) 
      stop("Must have initial.values$alpha.b > 0")
    
    if (initial.values$alpha.b < priors$alpha[1] | initial.values$alpha.b > priors$alpha[2])
      stop("invalid value for initial.values$alpha.b; must be between first and second elements of priors$alpha")
  }
  
  ##### remaining initializing values #####
  
  ###### delta ######
  
  # check presence
  if (!("delta" %in% names(initial.values))) {
    stop("initial.values must have element named delta")
  }
  
  # check format
  if (length(initial.values$delta) == 1) {
    initial.values$delta <- rep(initial.values$delta, lowertcount)
  }
  
  if (length(initial.values$delta) != lowertcount | !is.null(dim(initial.values$delta))) 
    stop("initial.values$delta must be a vector of length ncol(y) * (ncol(y) - 1)/2")
  
  ###### tau ######
  
  # only required if algorithm = SMN
  
  if (algorithm == "SMN") {
    
    # check presence
    if (!("tau" %in% names(initial.values))) {
      stop("When algorithm = SMN, initial.values must have element named tau")
    }
    
    # check format
    if (length(initial.values$tau) == 1) 
      initial.values$tau <- rep(initial.values$tau, lowertcount)
    
    if (length(initial.values$tau) != lowertcount) 
      stop("initial.values$tau must be a vector of length ncol(y) * (ncol(y) - 1)/2")
    
    # check validity
    if (!all(initial.values$tau > 0)) 
      stop("all values in initial.values$tau must be greater than zero")
  }
  
  ###### alpha.d ######
  
  # check presence
  if (!("alpha.d" %in% names(initial.values))) {
    stop("initial.values must have element named alpha.d")
  }
  
  if (length(initial.values$alpha.d) < 1) {
    stop("initial.values$alpha.d must be a scalar")
  } else if (length(initial.values$alpha.d) > 1) {
    initial.values$alpha.d <- initial.values$alpha.d[1] # ignore all but first element if more than value provided
  }
  
  # validity
  if (initial.values$alpha.d <= 0) 
    stop("Must have initial.values$alpha.d > 0")
  
  if (initial.values$alpha.d < priors$alpha[1] | initial.values$alpha.d > priors$alpha[2])
    stop("invalid value for initial.values$alpha.d; must be between first and second elements of priors$alpha")
  
  ###### gamma ######
  
  # check presence
  if (!("gamma" %in% names(initial.values))) {
    stop("initial.values must have element named gamma")
  }
  
  if (length(initial.values$gamma) != q | !is.null(dim(initial.values$gamma))) 
    stop("initial.values$gamma must be a vector of length ncol(y)")
  if (any(initial.values$gamma <= 0)) 
    stop("initial.values$gamma values must be greater than 0")
  
  #### addressing gamma prior ####
  
  # the only prior that does not have to be user-specified
  
  # check presence
  if (!("gamma" %in% names(priors))) {
    priors$gamma <- gamma.mom(initial.values$gamma)
  }
  
  # check format
  if (!is.null(dim(priors$gamma)))
    stop("priors$gamma should be a numeric vector, not a matrix")
  
  if (length(priors$gamma) != 2) 
    stop("priors$gamma should be vector of length 2")
  
  # check validity
  if (any(priors$gamma <= 0)) 
    stop("All values in priors$gamma must be greater than zero")
  
  #### preparation ####
  pos <- lapply(2:q, function(c) pos.vechrow(c, 1:(c-1))) # indexing for each delta_c, list of length q-1
  iter.1 <- iter + 1 # will include initial values as "step 0"
  
  # create .Random.seed if it does not exist
  if (!exists(".Random.seed")) {
    rnorm(1)
  }
  starting.seed <- .Random.seed
  
  if (algorithm == "MH") {
    
    if (debug) {
      if (meanzero) {
        out <- gmcb_mh_meanzero_debug(y = y, 
                                      d_init = initial.values$delta,
                                      gamma_init = initial.values$gamma, 
                                      alpha_d_init = initial.values$alpha.d, 
                                      gamma_prior = priors$gamma, 
                                      alpha_prior = priors$alpha, 
                                      tau_prior = priors$tau, 
                                      iter = iter, 
                                      d_scale = rw.scale$delta,
                                      alpha_d_scale = rw.scale$alpha.d, 
                                      pos = pos)
      } else {
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
      }
    } else {
      if (meanzero) {
        out <- gmcb_mh_meanzero_nodebug(y = y, 
                                        d_init = initial.values$delta,
                                        gamma_init = initial.values$gamma, 
                                        alpha_d_init = initial.values$alpha.d, 
                                        gamma_prior = priors$gamma, 
                                        alpha_prior = priors$alpha, 
                                        tau_prior = priors$tau, 
                                        iter = iter, 
                                        d_scale = rw.scale$delta,
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
    }
    
  } else {
    
    if (debug) {
      if (meanzero) {
        out <- gmcb_smn_meanzero_debug(y = y,
                                       d_init = initial.values$delta, 
                                       gamma_init = initial.values$gamma, 
                                       alpha_d_init = initial.values$alpha.d, 
                                       tau_init = initial.values$tau,
                                       gamma_prior = priors$gamma, 
                                       alpha_prior = priors$alpha, 
                                       tau_prior = priors$tau, 
                                       iter = iter, 
                                       alpha_d_scale = rw.scale$alpha.d, 
                                       pos = pos)
      } else {
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
      }
      
    } else {
      if(meanzero) {
        out <- gmcb_smn_meanzero_nodebug(y = y,
                                         d_init = initial.values$delta, 
                                         gamma_init = initial.values$gamma, 
                                         alpha_d_init = initial.values$alpha.d, 
                                         tau_init = initial.values$tau,
                                         gamma_prior = priors$gamma, 
                                         alpha_prior = priors$alpha, 
                                         tau_prior = priors$tau, 
                                         iter = iter, 
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
      seconds.per.step[t] <- seconds.per.step[t]/86400 # convert to days
    }
  }
  timing <- paste0(round(seconds.per.step, digits = 4), units.per.step)
  names(timing) <- c("Pre-sampling computations", "Sampling")
  
  ## format starting values for continuing runs
  if (meanzero) {
    final.values <- list(delta = out$mcmc$delta[iter.1,],
                         gamma = out$mcmc$gamma[iter.1,], 
                         alpha.d = out$mcmc$alpha_d[iter.1])
    
    if (algorithm == "SMN") {
      final.values <- c(final.values, list(tau = out$mcmc$tau[iter.1,]))
    } 
  } else {
    final.values <- list(b = out$mcmc$b[iter.1,], delta = out$mcmc$delta[iter.1,],
                         gamma = out$mcmc$gamma[iter.1,], alpha.b = out$mcmc$alpha_b[iter.1],
                         alpha.d = out$mcmc$alpha_d[iter.1])
    
    if (algorithm == "SMN") {
      final.values <- c(final.values, list(lambda = out$mcmc$lambda[iter.1,],
                                           tau = out$mcmc$tau[iter.1,]))
    } 
  }
  
  final.out <- list(mcmc = out$mcmc, acceptances = out$acceptances,
                    start.seed = starting.seed, end.seed = ending.seed,
                    timing = timing, final.values = final.values,
                    y = y, x = x, meanzero = meanzero,
                    initial.values = initial.values,
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

gmcb.gbridge <- function(y, x = NULL, meanzero = FALSE,
                         initial.values = NULL,
                         priors,
                         rw.scale,
                         iter = 1e4, 
                         algorithm = c("MH", "SMN"),
                         debug = FALSE) {
  
  if (missing(meanzero)) meanzero <- y$meanzero
  if (missing(priors)) priors <- y$priors
  if (missing(rw.scale)) rw.scale <- y$rw.scale
  if (missing(iter)) iter <- y$iter
  if (missing(algorithm)) algorithm <- y$algorithm
  if (missing(debug)) debug <- y$debug
  
  response <- y$y
  x <- y$x
  
  initial.values <- y$final.values # final values from previous run are now initializing values
  
  assign(".Random.seed", y$end.seed, .GlobalEnv)
  
  out <- gmcb.matrix(response, x = x, meanzero = meanzero,
                     initial.values = initial.values,
                     priors = priors,
                     rw.scale = rw.scale,
                     iter = iter, 
                     algorithm = algorithm,
                     debug = debug)
  
  # since initialized using final values from previous run, remove initial values from matrix output
  if (meanzero) {
    out$mcmc$delta <- out$mcmc$delta[-1,,drop = FALSE]
    out$mcmc$gamma <- out$mcmc$gamma[-1,,drop = FALSE]
    out$mcmc$tau <- out$mcmc$tau[-1,,drop = FALSE]
    out$mcmc$alpha_d <- out$mcmc$alpha_d[-1]
  } else {
    out$mcmc$b <- out$mcmc$b[-1,,drop = FALSE]
    out$mcmc$delta <- out$mcmc$delta[-1,,drop = FALSE]
    out$mcmc$gamma <- out$mcmc$gamma[-1,,drop = FALSE]
    out$mcmc$lambda <- out$mcmc$lambda[-1,,drop = FALSE]
    out$mcmc$tau <- out$mcmc$tau[-1,,drop = FALSE]
    out$mcmc$alpha_b <- out$mcmc$alpha_b[-1]
    out$mcmc$alpha_d <- out$mcmc$alpha_d[-1]
  }
  
  return(out)
}
  