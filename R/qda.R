
qda <- function(test, classes, class.prior, class.means, class.precisions,
                parallel = 0) {
  # test: matrix of test observations, m x q
  # classes: class labels, eg. 0, 1
  # class.prior: the prior probability of training observations falling in each class
  # class.means: list of class means, with each mean a vector of length q
  # class.precisions: list of precision matrices for each class
  # parallel: nonnegative integer value. If 0 or 1, no parallelization, otherwise 
    # if the value is larger, will assume that is the number of cores to use for parallelization
    # if not an integer value, will take the floor of the value
  
  # will assume that the means and precisions provided are in order of the classes 
  # specified in classes
  
  #### error checking ####
  stopifnot(length(classes) == length(class.prior))
  stopifnot(length(classes) == length(class.means))
  stopifnot(length(classes) == length(class.precisions))
  stopifnot(is.numeric(parallel))
  stopifnot(parallel >= 0)
  parallel <- floor(parallel)
  
  if (parallel < 2) {
    test.class <- numeric(nrow(test))
    for (i in 1:nrow(test)) {
      d.f <- numeric(length(classes))
      for (j in 1:length(classes)) {
        d.f[j] <- determinant(class.precisions[[j]], logarithm = TRUE)$modulus/2 -
          (1/2)* mahalanobis(x = test[i,], center = class.means[[j]], cov = class.precisions[[j]],
                             inverted = TRUE) + log(class.prior[j])
      }
      test.class[i] <- classes[which.max(d.f)]
    }
  } else {
    clus <- makeCluster(parallel)
    registerDoParallel(clus)
    
    test.class <- foreach(i = 1:nrow(test), .combine = c) %dopar% {
      d.f <- numeric(length(classes))
      for (j in 1:length(classes)) {
        d.f[j] <- determinant(class.precisions[[j]], logarithm = TRUE)$modulus/2 -
          (1/2)* mahalanobis(x = test[i,], center = class.means[[j]], cov = class.precisions[[j]],
                             inverted = TRUE) + log(class.prior[j])
      }
      return(classes[which.max(d.f)])
    }
  }
  
  stopCluster(clus)
  
  test.class
}

