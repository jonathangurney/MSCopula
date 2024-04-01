#' ptCopula
#'
#' A function to compute the distribution function of the Student-t copula.
#'
#' This function computes the Student-t copula distribution function at \eqn{n}
#' points, where the \eqn{i}th point in \eqn{d} dimensions is specified by the
#' \eqn{i}th row of the matrix \code{U}. \code{rho} specifies the correlation
#' matrix used to parameterise the copula distribution function.
#'
#' @param U an \eqn{n x d} matrix giving \eqn{n} points in \eqn{[0, 1]^d} at
#' which to evaluate the distribution function.
#' @param rho a \eqn{d x d} correlation matrix, a vector of length \eqn{n} or a
#' list of \eqn{n} \eqn{d x d} correlation matrices. \code{rho} can be a vector
#' of scalar values only if \eqn{d = 2}.
#' @param df a scalar or a numeric vector of length \eqn{n} determining the degrees
#' of freedom of the Student-t distribution used in the copula.
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
#' df <- 3
#' rho <- c(0.4, 0.1, 0.7)
#' U <- matrix(runif(n * d), ncol = d)
#' cdft <- ptCopula(U, rho, df)
ptCopula <- function(U, rho, df) {
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

  # In the case that rho is a vector of scalar values, extend it if the length
  # of rho is < n and trim it is the length of rho is > n and coerce to matrix form
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

  # Create vector of degrees of freedom
  df <- rep_len(df, length.out = n)

  Z <- qt(U, df = df)

  cdf_eval <- function(i, Z, Rho, df) {
    return(mvtnorm::pmvt(lower = -Inf, upper = Z[i, ], corr = Rho[[i]], df = df[i])[1])
  }

  Y <- sapply(1:n, FUN = cdf_eval, Z = Z, Rho = Rho, df = df)
  return(Y)
}

#' dtCopula
#'
#' A function to compute the density function of the Student-t copula.
#'
#' This function computes the Student-t copula density at \eqn{n} points, where
#' the \eqn{i}th point in \eqn{d} dimensions is specified by the \eqn{i}th row
#' of the matrix \code{U}. \code{rho} specifies the correlation matrices used to
#' parameterise the copula density function.
#'
#' @param U an \eqn{n x d} matrix giving \eqn{n} points in \eqn{[0, 1]^d} at
#' which to evaluate the distribution function.
#' @param rho a \eqn{d x d} correlation matrix, a vector of length \eqn{n} or a
#' list of \eqn{n} \eqn{d x d} correlation matrices. \code{rho} can be a vector
#' of scalar values only if \eqn{d = 2}.
#' @param df a scalar or a numeric vector of length \eqn{n} determining the degrees
#' of freedom of the Student-t distribution used in the copula.
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
#' df <- 3
#' rho <- replicate(n, p2rho(runif(d*(d - 1)/2, max = 0.6)), simplify = FALSE)
#' U <- matrix(runif(n * d), ncol = d)
#' pdft <- dtCopula(U, rho, df)
dtCopula <- function(U, rho, df) {
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

  # In the case that rho is a vector of scalar values, extend it if the length
  # of rho is < n and trim it is the length of rho is > n and coerce to matrix form
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

  # Create vector of degrees of freedom
  df <- rep_len(df, length.out = n)

  Z <- qt(U, df = df)

  dens_fun <- function(x, rho, df) {
    return(mvtnorm::dmvt(x, sigma = rho, df = df, log = FALSE)/prod(dt(x, df = df, log = FALSE)))
  }

  dens_fun_eval <- function(i, Z, Rho, df) {
    return(dens_fun(Z[i, ], Rho[[i]], df[i]))
  }

  Y <- sapply(1:n, FUN = dens_fun_eval, Z = Z, Rho = Rho, df = df)
  return(Y)
}

#' tLogLik
#'
#' A function to compute the log-likelihood function for parameters of the
#' Student-t copula.
#'
#' This function computes the Student-t copula log-likelihood at \eqn{n} points, where
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
#' @param df a scalar or a numeric vector of length \eqn{n} determining the degrees
#' of freedom of the Student-t distribution used in the copula.
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
#' df <- 3
#' rho <- replicate(n, p2rho(runif(d*(d - 1)/2, max = 0.6)), simplify = FALSE)
#' U <- matrix(runif(n * d), ncol = d)
#' loglikGaussian <- GaussianLogLik(rho, U, df)
tLogLik <- function(rho, U, df) {
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

  # In the case that rho is a vector of scalar values, extend it if the length
  # of rho is < n and trim it is the length of rho is > n and coerce to matrix form
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

  loglik <- log(dtCopula(U, rho, df))
  loglik[is.infinite(loglik)] <- -1e9
  return(loglik)
}

