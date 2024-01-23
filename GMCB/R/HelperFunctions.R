
#### Convenience functions ####

# a function that converts a q x q matrix to a vector by stacking the
# lower-triangular elements by row
vech.row <- function(mat, keep.diag = FALSE) {
  # diag = FALSE does not keep the diagonals and returns a vector with length q(q-1)/2 
  # diag = TRUE keeps the diagonals and returns a vector with length q(q+1)/2
  temp <- mat
  temp[upper.tri(temp, diag = !keep.diag)] <- NA
  return(as.numeric(na.omit(matrixcalc::vec(t(temp)))))
}

# computing position of lower-triangular elements in row vech
pos.vechrow <- function(c, j) {
  if (c == 2) {
    return(1)
  } else {
    return(sum(2:(c-1)-1) + j)
  }
}

#### functions for working with GMCB output ####

# function to convert an m x q(q-1)/2 matrix to a list of length m, with elements
# q x q unit lower-triangular matrices
delta.to.matrix <- function(delta, q) {
  # delta: m x q(q-1)/2 matrix
  # elements fill the lower-triangular matrix by row; that is, 
  # the first element of delta is element [2,1]
  # the second element is element [3,1]
  # the third element is element [3,2]
  # the fourth element is element [4,1]
  # etc
  # q: dimension of unit lower-triangular matrix
  m <- nrow(delta)
  out <- vector("list", length = m)
  for (i in 1:m) {
    out[[i]] <- delta_to_matrix_inner_smn(delta[i,], q)
  }
  return(out)
}

# function to convert an m x pq matrix to a list of length m with elements being p x q matrices
b.to.matrix <- function(b, p, q) {
  # b: m x pq matrix 
  # p: number of rows for matrices
  # q: number of columns for matrices
  m <- nrow(b)
  out <- vector("list", length = m)
  for (i in 1:m) {
    out[[i]] <- matrix(b[i,], p, q)
  }
  return(out)
}

# function to calculate the posterior samples for the covariance/precision matrix
# either list of posterior samples for the precision matrix
# or a list of two lists, one containing samples of the precision matrix, 
# the other the samples of the covariance matrix
prec.mat <- function(delta, gamma, cov = TRUE) {
  # delta: a m x q(q-1)/2 matrix, where m is the number of samples
  # gamma: a m x q matrix, where m is the number of samples
  # cov: argument stating whether the covariance matrix samples are also desired
  
  q <- ncol(gamma)
  stopifnot(ncol(delta) == q*(q-1)/2)
  
  delta.list <- delta.to.matrix(delta, q)
  prec.mat <- lapply(1:length(delta.list), 
                     function(x) t(delta.list[[x]]) %*% diag(1/gamma[x,]) %*% delta.list[[x]])
  if (cov) {
    inv.delta.list <- lapply(delta.list, solve, tol = 1e-30)
    cov.mat <- lapply(1:length(delta.list), 
                      function(x) inv.delta.list[[x]] %*% diag(gamma[x,]) %*% t(inv.delta.list[[x]]))
    return(list(prec.mat = prec.mat, cov.mat = cov.mat))
  } else {
    return(prec.mat)
  }
}

# function for turning a list of q x q matrices with length m to a m x q(q+1)/2 matrix
vech.list <- function(list) {
  # list: length assumed to be the number of observations, dimension of each element assumed to be the same
  
  m <- length(list)
  q <- nrow(list[[1]])
  vechs <- sapply(list, matrixcalc::vech)
  return(t(vechs))
}

#### General helper functions ####

# function for generating CS, AR(1), and independent covariance matrices and their inverses
matgen <- function(p, rho, sigma = 1, type = c("AR", "CS", "cI")) {
  # p = dimension
  # rho = correlation 
  # ignored if type = cI
  # sigma = common variance if AR or CS, diagonal values for cI
  # if sigma is a vector and type is AR or CS, only the first element is used
  # type = AR, CS, cI
  stopifnot(sigma > 0)
  stopifnot(p > 1)
  type <- match.arg(type)
  
  if (type == "AR") {
    stopifnot(abs(rho) <= 1)
    cov.mat <- matrix(nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in 1:p) {
        cov.mat[i,j] <- rho^abs(i - j)
      }
    }
    cov.mat <- cov.mat*sigma[1]
    
    # true precision matrix
    prec.mat <- matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in 1:p) {
        if (i == j) {
          if (i %in% c(1, p)) {
            prec.mat[i,j] <- 1/(1 - rho^2) # the (1,1) and (p,p) elements are 1/(1-rho^2)
          } else {
            prec.mat[i,j] <- (1 + rho^2)/(1 - rho^2)
          }
        } else if ((j == i - 1) | (j == i + 1)) {
          prec.mat[i,j] <- -rho/(1 - rho^2)
        }
      }
    }
    prec.mat <- 1/sigma[1] * prec.mat
  } else if (type == "CS") { 
    stopifnot(abs(rho) <= 1)
    cov.mat <- matrix(rho, nrow = p, ncol = p)
    diag(cov.mat) <- 1
    cov.mat <- sigma[1]*cov.mat
    
    prec.mat <- solve(cov.mat)
  } else {
    if (length(sigma) == 1)
      sigma <- rep(sigma, length.out = p)
    cov.mat <- diag(sigma)
    prec.mat <- diag(1/sigma)
  }
  return(list(cov.mat = cov.mat, prec.mat = prec.mat))
}

# function for computing the matrices T and D
mcd <- function(cov.mat) {
  # cov.mat: the covariance matrix corresponding to the precision matrix we are 
  # finding the modified Cholesky decomposition for
  
  q <- nrow(cov.mat)
  l <- t(chol(cov.mat))
  d <- matrix(0, nrow = q, ncol = q)
  diag(d) <- diag(l)^2
  
  cinv <- diag(1/diag(l))
  tinv <- l %*% cinv
  tt <- forwardsolve(tinv, diag(q))
  
  return(list(t = tt, d = d))
}