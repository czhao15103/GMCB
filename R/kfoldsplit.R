
#### Splitting for K-fold CV ####

kfoldsplit <- function(n, folds = 5) {
  # n = number of observations
  # folds = number of disjoint folds desired
  
  samples <- 1:n
  
  leftover <- n %% folds 
  
  obsinfolds <- vector("list", folds)
  if (leftover == 0) {
    foldn <- n/folds # set number of folds per observation
    for (i in 1:folds) {
      obsinfolds[[i]] <- sample(samples, size = foldn)
      samples <- setdiff(samples, obsinfolds[[i]])
    }
  } else {
    foldn <- rep(floor(n/folds), folds) # vector of fold sizes when n not divisible by folds
    foldn[1:leftover] <- foldn[1:leftover] + 1
    for (i in 1:folds) {
      obsinfolds[[i]] <- sample(samples, size = foldn[i])
      samples <- setdiff(samples, obsinfolds[[i]])
    }
  }
  return(obsinfolds)
}
