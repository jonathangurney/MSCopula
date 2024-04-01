#' pGaussianCopula
#'
#' A function to compute the distribution function of the Gaussian copula.
#'
#' This function computes the Gaussian copula distribution function at \eqn{n}
#' points, where the \eqn{i}th point in \eqn{d} dimensions is specified by the
#' \eqn{i}th row of the matrix \code{U}. \code{rho} specifies the correlation
#' matrix used to parameterise the copula distribution function.
#'
#' @param U an \eqn{n x d} matrix giving \eqn{n} points in \eqn{[0, 1]^d} at
#' which to evaluate the distribution function.
#' @param rho either a \eqn{d x d} symmetric, positive definite matrix of
#' correlation coefficients or a vector of length \eqn{d*(d-1)/2}. If \code{rho}
#' is a vector, the correlation matrix is constructed by filling the upper
#' triangle of the matrix row-by-row with elements of \code{rho}.
#'
#' @import matrixcalc mvtnorm
#'
#' @return The output is an \eqn{n x 1} vector of cdf values.
#'
#' @references
#' \insertref{nelsen2006}{MSCopula}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' d <- 3
#' rho <- c(0.4, 0.1, 0.7)
#' U <- matrix(runif(n * d), ncol = d)
#' cdfGaussian <- pGaussianCopula(U, rho)
pGaussianCopula <- function(U, rho) {
  # Check if U is a matrix and coerce to matrix form if necessary. We assume that
  # U is a (1 x d) vector in this instance.
  if (!any(class(U) %in% c('matrix'))) {
    U <- matrix(U, nrow = 1)
  }

  n <- nrow(U)
  d <- ncol(U)

  if (d < 2) {
    stop('Copula must have at least 2 dimensions.')
  }

  # If rho has been given as a vector of coefficients coerce to matrix form
  if (!(any(class(rho) %in% c("matrix")))) {
    if (length(rho) != (d * (d - 1)/2)) {
      stop('Number of correlation coefficients does not match dimensions of U')
    }

    rho <- p2rho(rho)
  }

  # Check that rho is a symmetric, positive definite matrix with elements in
  # [-1, 1]
  if (!(matrixcalc::is.symmetric.matrix(rho))) {
    stop('Correlation matrix not symmetric')
  }
  if (!(matrixcalc::is.positive.definite(rho))) {
    stop('Correlation matrix not positive definite')
  }

  if (any(abs(rho) > 1)) {
    stop('Correlation matrix contains element outside interval [-1, 1]')
  }

  X <- qnorm(U)
  Y <- apply(X, FUN = mvtnorm::pmvnorm, MARGIN = 1, corr = rho, lower = -Inf)
  return(Y)
}

#' dGaussianCopula
#'
#' A function to compute the density function of the Gaussian copula.
#'
#' This function computes the Gaussian copula density at \eqn{n} points, where
#' the \eqn{i}th point in \eqn{d} dimensions is specified by the \eqn{i}th row
#' of the matrix \code{U}. \code{rho} specifies the correlation matrices used to
#' parameterise the copula density function.
#'
#' @param U an \eqn{n x d} matrix giving \eqn{n} points in \eqn{[0, 1]^d} at
#' which to evaluate the distribution function.
#' @param rho a \eqn{d x d} correlation matrix, a vector of length \eqn{n} or a
#' list of \eqn{n} \eqn{d x d} correlation matrices. \code{rho} can be a vector
#' of scalar values only if \eqn{d = 2}.
#'
#' @import matrixcalc mvtnorm
#'
#' @return The output is an \eqn{n x 1} vector of cdf values.
#'
#' @references
#' \insertref{nelsen2006}{MSCopula}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' d <- 3
#' rho <- replicate(n, p2rho(runif(d*(d - 1)/2, max = 0.6)), simplify = FALSE)
#' U <- matrix(runif(n * d), ncol = d)
#' pdfGaussian <- dGaussianCopula(U, rho)
dGaussianCopula <- function(U, rho) {
  # Check if U is a matrix and coerce to matrix form if necessary. We assume that
  # U is a (1 x d) vector in this instance.
  if (!any(class(U) %in% c('matrix'))) {
    U <- matrix(U, nrow = 1)
  }

  n <- nrow(U)
  d <- ncol(U)

  if (d < 2) {
    stop('Copula must have at least 2 dimensions')
  }

  # In the case that rho is a vector of scalar values, extend it if the length
  # of rho is < n and trim it is the length of rho is > n and coerce
  if (any(class(rho) %in% c('numeric'))) {
    if (!(d == 2)) {
      stop('A vector of correlation coefficients can only be given if the copula is bivariate')
    }

    # Make Rho the correct length
    if (length(rho) < n) {
      message('rho is of length < n and has been extended by recycling values.')
      rho <- rep_len(rho, length.out = n)
    } else if (length(rho) > n) {
      message('rho is of length > n and the first n items have been used.')
      rho <- rho[1:n]
    }

    # Coerce rho to be a list of correlation matrices
    Rho <- lapply(rho, p2rho)
  }

  # If rho is a matrix, create a list of n replications of rho
  if (any(class(rho) %in% c('matrix'))) {
    Rho <- replicate(n, rho, simplify = FALSE)
  }

  # If a list of matrices has been supplied set Rho equal to rho,
  if (any(class(rho) %in% c('list'))) {
    if (length(rho) < n) {
      message('rho has length < n and has been extended by recycling values.')
      Rho <- rep_len(rho, length.out = n)
    } else if (length(rho) > n) {
      message('rho is of length > n and the first n items have been used.')
      Rho <- rho[1:n]
    } else {
      Rho <- rho
    }
  }

  # Check that all correlation matrices are positive definite
  is_pos_def <- sapply(Rho, matrixcalc::is.positive.definite)
  if (!all(is_pos_def)) {
    stop('Not all correlation matrices supplied are positive definite.')
  }

  # Check that all correlation matrices are symmetric
  is_symm <- sapply(Rho, matrixcalc::is.symmetric.matrix)
  if (!all(is_symm)) {
    stop('Not all correlation matrices supplied are symmetric.')
  }

  dens_fun <- function(x, rho) {
    return(mvtnorm::dmvnorm(x, sigma = rho)/prod(dnorm(x)))
  }

  dens_fun_index <- function(i) {
    return(dens_fun(X[i, ], Rho[[i]]))
  }

  X <- qnorm(U)
  Y <- sapply(1:n, FUN = dens_fun_index)
  return(Y)
}

#' GaussianLogLik
#'
#' A function to compute the log-likelihood function for parameters of the
#' Gaussian copula.
#'
#' This function computes the Gaussian copula log-likelihood at \eqn{n} points, where
#' the \eqn{i}th point in \eqn{d} dimensions is specified by the \eqn{i}th row
#' of the matrix \code{U}. \code{rho} specifies the correlation matrices used to
#' parameterise the copula density function.
#'
#' @param rho a \eqn{d x d} correlation matrix, a vector of length \eqn{n} or a
#' list of \eqn{n} \eqn{d x d} correlation matrices. \code{rho} can be a vector
#' of scalar values only if \eqn{d = 2}. If \code{rho} has length less than \eqn{n}
#' then the values in \code{rho} are recycled to give \code{rho} length \eqn{n}
#' @param U an \eqn{n x d} matrix giving \eqn{n} points in \eqn{[0, 1]^d} at
#' which to evaluate the log-likelihood.
#'
#' @return The output is an \eqn{n x 1} vector of log-likelihood values. Any
#' value of the log-likelihood which evaluates to infinity is replaced with -1e9.
#'
#' @references
#' \insertref{nelsen2006}{MSCopula}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' d <- 3
#' rho <- replicate(n, p2rho(runif(d*(d - 1)/2, max = 0.6)), simplify = FALSE)
#' U <- matrix(runif(n * d), ncol = d)
#' loglikGaussian <- GaussianLogLik(rho, U)
GaussianLogLik <- function(rho, U) {
  # Check if U is a matrix and coerce to matrix form if necessary. We assume that
  # U is a (1 x d) vector in this instance.
  if (!any(class(U) %in% c('matrix'))) {
    U <- matrix(U, nrow = 1)
  }

  loglik <- log(dGaussianCopula(U, rho))
  loglik[is.infinite(loglik)] <- -1e9
  return(loglik)
}

# A helper function that converts a numeric vector to a correlation matrix.
p2rho <- function(p) {
  d <- 1/2 + sqrt(1/4 + 2*length(p))
  if (d %% 1 != 0) {
    stop('p is not of the correct length.')
  }
  rho <- matrix(rep(0, d*d), nrow = d, ncol = d)
  rho[upper.tri(rho)] <- p
  rho <- t(rho) + rho
  diag(rho) <- 1
  return(rho)
}

# Does this function even get used anywhere?
CorrMat <- function(n, d) {
  fun <- function(i, d) {
    m <- 0.5*d*(d - 1)
    cor.mat <- diag(d)
    cor.mat[upper.tri(cor.mat)] <- vals[((i - 1)*m + 1):(i*m)]
    cor.mat[lower.tri(cor.mat)] <- vals[((i - 1)*m + 1):(i*m)]
    return(cor.mat)
  }

  vals <- runif(n*0.5*d*(d-1), max = 0.725)
  mat_list <- lapply(X = 1:n, FUN = fun, d = d)
  return(mat_list)
}
