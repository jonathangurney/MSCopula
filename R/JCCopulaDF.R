#' pJCCopula
#'
#' A function to compute the distribution function of the Joe-Clayton copula.
#'
#' This function computes the Joe-Clayton copula distribution function at \eqn{n}
#' points, where the \eqn{i}th point in \eqn{d} dimensions is specified by the
#' \eqn{i}th row of the matrix \code{U}. If \code{theta} is given as a numeric
#' vector then these parameter values are used to compute the cdf at all points.
#' If \code{theta} is a matrix, then the parameters in the \eqn{i}th row of
#' \code{theta} are used to compute the cdf at the point specified in the \eqn{i}th
#' row of \code{U}.
#'
#' @param U an \eqn{n x d} matrix giving n points in \eqn{[0, 1]^d} at which to evaluate
#' the cdf of the Joe-Clayton copula.
#' @param theta either a numeric vector of length 2 or an \eqn{n x 2} matrix of
#' parameter values in \eqn{[0, 1]}. Values in the first column of \code{theta}
#' determine the lower tail dependence, while those in the second column determine
#' the upper tail dependence. If \code{theta} is given as a matrix, values in the
#' \eqn{i}th row of \code{theta} are used to compute the cdf at the point specified
#' by row \eqn{i} of \code{U}. If \code{theta} is given as a numeric vector of
#' length 2 then these values are used to parameterise the cdf at all points.
#'
#' @return The output is an \eqn{n x 1} vector of cdf values.
#'
#' @export
#'
#' @references
#' \insertRef{nelsen2006}{MSCopula}
#' \insertRef{patton2006}{MSCopula}
#'
#' @importFrom Rdpack reprompt
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' d <- 3
#' theta <- cbind(runif(n), runif(n))
#' U <- matrix(runif(n * d), ncol = d)
#' cdfJoeClayton <- pJCCopula(U, theta)
pJCCopula <- function(U, theta) {
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


  if (any(class(theta) %in% c('numeric'))) {
    # Check that the correct number of parameters has been specified
    if (length(theta) != 2) {
      stop('Two parameters required to compute distribution function.')
    }

    # Check that the parameters lie in [0, 1]
    if ((any(theta < 0)) | (any(theta > 1))) {
      stop('Parameters must lie in the interval [0, 1].')
    }

    # Convert theta to a matrix with n rows
    theta <- matrix(rep(theta, each = n), ncol = 2)
  } else if (any(class(theta) %in% c('matrix'))) {
    # Check that the matrix is of the correct dimensions
    if (ncol(theta) != 2) {
      stop('Parameter matrix theta does not have correct dimensions.')
    }

    # If theta has fewer than n rows recycle values so that theta has exactly n
    # rows
    if (nrow(theta) < n) {
      theta1 <- rep_len(theta[, 1], length.out = n)
      theta2 <- rep_len(theta[, 2], length.out = n)

      theta <- cbind(theta1, theta2)
    }

    # If theta has more than n rows choose only the first n
    if (nrow(theta) > n) {
      theta <- theta[1:n, ]
    }

    # Check that all parameter values are in correct range
    if ((any(theta < 0)) | (any(theta > 1))) {
      stop('Parameter values must lie in the interval [0, 1].')
    }
  }

  # Define function that applies computes the distribution function at each point in U
  cdf_eval_fun <- function(i, U, theta) {
    return(JCCDFeval(U[i, ], theta[i, ]))
  }

  # Apply cdf_eval_fun over all rows of U and theta
  Y <- sapply(X = 1:n, FUN = cdf_eval_fun, U = U, theta = theta)
  return(Y)
}

#' dJCCopula
#'
#' A function to compute the density function of the Joe-Clayton copula.
#'
#' This function computes the Joe-Clayton copula density function at \eqn{n}
#' points, where the \eqn{i}th point in \eqn{d} dimensions is specified by the
#' \eqn{i}th row of the matrix \code{U}. If \code{theta} is given as a numeric
#' vector then these parameter values are used to compute the density at all points.
#' If \code{theta} is a matrix, then the parameters in the \eqn{i}th row of
#' \code{theta} are used to compute the pdf at the point specified in the \eqn{i}th
#' row of \code{U}.
#'
#' @param U an \eqn{n x d} matrix giving n points in \eqn{[0, 1]^d} at which to evaluate
#' the density of the Joe-Clayton copula.
#' @param theta either a numeric vector of length 2 or an \eqn{n x 2} matrix of
#' parameter values in \eqn{[0, 1]}. Values in the first column of \code{theta}
#' determine the lower tail dependence, while those in the second column determine
#' the upper tail dependence. If \code{theta} is given as a matrix, values in the
#' \eqn{i}th row of \code{theta} are used to compute the pdf at the point specified
#' by row \eqn{i} of \code{U}. If \code{theta} is given as a numeric vector of
#' length 2 then these values are used to parameterise the density at all points.
#'
#' @return The output is an \eqn{n x 1} vector of pdf values. If any point returns
#' a value of the density that evaluates to \code{Inf} then NaN is returned.
#'
#' @references
#' \insertRef{nelsen2006}{MSCopula}
#' \insertRef{patton2006}{MSCopula}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' d <- 3
#' theta <- cbind(runif(n), runif(n))
#' U <- matrix(runif(n * d), ncol = d)
#' pdfJoeClayton <- dJCCopula(U, theta)
dJCCopula <- function(U, theta) {
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


  if (any(class(theta) %in% c('numeric'))) {
    # Check that the correct number of parameters has been specified
    if (length(theta) != 2) {
      stop('Two parameters required to compute distribution function.')
    }

    # Check that the first parameter lies in the interval [0, \infty)
    if ((any(theta < 0)) | (any(theta > 1))) {
      stop('Parameter values outside of correct range.')
    }

    # Convert theta to a matrix with n rows
    theta <- matrix(rep(theta, each = n), ncol = 2)
  } else if (any(class(theta) %in% c('matrix'))) {
    # Check that the matrix is of the correct dimensions
    if (ncol(theta) != 2) {
      stop('Parameter matrix theta does not have correct dimensions.')
    }

    # If theta has fewer than n rows recycle values so that theta has exactly n
    # rows
    if (nrow(theta) < n) {
      theta1 <- rep_len(theta[, 1], length.out = n)
      theta2 <- rep_len(theta[, 2], length.out = n)

      theta <- cbind(theta1, theta2)
    }

    # If theta has more than n rows choose only the first n
    if (nrow(theta) > n) {
      theta <- theta[1:n, ]
    }

    # Check that all parameter values are in correct range
    if ((any(theta < 0)) | (any(theta > 1))) {
      stop('Parameter values not in correct range.')
    }
  }

  pdf_eval_fun <- function(i, U, theta, dim) {
    return(JCPDFeval(U[i, ], dim, theta[i, ]))
  }

  Y <- sapply(1:n, FUN = pdf_eval_fun, U = U, theta = theta, dim = d)
  Y[is.infinite(Y)] <- NaN
  return(Y)
}

#' JCLogLik
#'
#' A function to compute the log-likelihood function of the Joe-Clayton copula.
#'
#' This function computes the Joe-Clayton copula log-likelihood at \eqn{n}
#' points, where the \eqn{i}th point in \eqn{d} dimensions is specified by the
#' \eqn{i}th row of the matrix \code{U}. If \code{theta} is given as a numeric
#' vector then these parameter values are used to compute the density at all points.
#' If \code{theta} is a matrix, then the parameters in the \eqn{i}th row of
#' \code{theta} are used to compute the pdf at the point specified in the \eqn{i}th
#' row of \code{U}.
#'
#' @param theta either a numeric vector of length 2 or an \eqn{n x 2} matrix of
#' parameter values in \eqn{[0, 1]}. Values in the first column of \code{theta}
#' determine the lower tail dependence, while those in the second column determine
#' the upper tail dependence. If \code{theta} is given as a matrix, the \eqn{i}th row
#' specifies the point at which the log-likelihood is computed using the data in
#' row \eqn{i} of the matrix \code{U}.
#' @param U an \eqn{n\times d}{n x d} matrix giving n points in \eqn{[0, 1]^d} at which to evaluate
#' the log-likelihood of the Joe-Clayton copula.
#'
#' @return The output is an \eqn{n x 1} vector of log-likelihood values. If any
#' point returns a value of the log-likelihood that evaluates to \code{Inf} then
#' -1e9 is returned.
#'
#' @references
#' \insertRef{nelsen2006}{MSCopula}
#' \insertRef{patton2006}{MSCopula}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' d <- 3
#' theta <- cbind(runif(n), runif(n))
#' U <- matrix(runif(n * d), ncol = d)
#' LogLikJoeClayton <- JCLogLik(theta, U)
JCLogLik <- function(theta, U) {
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

  Y <- log(dJCCopula(U, theta))
  Y[is.infinite(Y)] <- -1e9
}

JCCDFexpr <- function(dim) {
  CDFexpr <- paste0('(1 - (1 - (1 - u', 1:dim, ')^kappa)^(-gamma))', collapse = ' + ')
  CDFexpr <- paste0('1 - (1 - (1 - (', CDFexpr, '))^(-1/gamma))^(1/kappa)')
  CDFexpr <- parse(text = CDFexpr)
  return(CDFexpr)
}

JCCDFeval <- function(u, theta) {
  gamma <- -1/log2(theta[1])
  kappa <- 1/log2(2 - theta[2])
  dim <- length(u)
  CDFexpr <- JCCDFexpr(dim)
  for (i in 1:d) {
    assign(paste0('u', i), u[i])
  }
  y <- eval(CDFexpr)
  return(y)
}

JCPDFexpr <- function(dim) {
  CDFexpr <- JCCDFexpr(dim)
  PDFexpr <- CDFexpr
  for (i in 1:dim) {
    PDFexpr <- D(PDFexpr, paste0("u", i))
  }
  return(PDFexpr)
}

JCPDFeval <- function(u, dim, theta) {
  gamma <- -1/log2(theta[1])
  kappa <- 1/log2(2 - theta[2])
  pdfexpr <- JCPDFexpr(dim)
  for (i in 1:dim) {
    assign(paste0("u", i), u[i])
  }
  y <- eval(pdfexpr)
  return(y)
}
