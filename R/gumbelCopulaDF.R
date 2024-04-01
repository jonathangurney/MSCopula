#' pGumbelCopula
#'
#' A function to compute the distribution function of the Gumbel copula.
#'
#' This function computes the Gumbel copula distribution function at \eqn{n}
#' points, where the \eqn{i}th point in \eqn{d} dimensions is specified by the
#' \eqn{i}th row of the matrix \code{U}. The copula dependence parameter used to
#' compute the distribution function at each point is the scalar \code{theta}.
#'
#' @param U an \eqn{n x d} matrix giving n points in \eqn{[0, 1]^d} at which to evaluate
#' the cdf of the Gumbel copula.
#' @param theta scalar dependence parameter for the Gumbel copula. Restricted
#' to the interval \eqn{[1, \infty)}.
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
#' theta <- runif(1, min = 0, max = 10)
#' U <- matrix(runif(n * d), ncol = d)
#' cdfGumbel <- pGumbelCopula(U, theta)
pGumbelCopula <- function(U, theta) {
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

  if (length(theta) > 1) {
    stop('Parameter is not a scalar')
  }

  if (theta < 1) {
    stop('Parameter not in valid range')
  }

  X <- (-log(U))^theta
  Y <- apply(X, FUN = sum, MARGIN = 1)
  Y <- exp(-Y^(1/theta))
  return(Y)
}

#' dGumbelCopula
#'
#' A function to compute the Gumbel copula density.
#'
#' This function computes the Gumbel copula density function at \eqn{n} points,
#' where the \eqn{i}th point in \eqn{d} dimensions is specified by the \eqn{i}th
#' row of the matrix \code{U}. The copula dependence parameter used to compute
#' the density at each point is the \eqn{i}th value in the vector \code{theta}.
#' If \code{theta} is of length \eqn{< n}, values from \code{theta} are recycled
#' to produce a vector of length \eqn{n}
#'
#' @param U an \eqn{n x d} matrix giving n points in \eqn{[0, 1]^d} at which to evaluate
#' the pdf of the Gumbel copula
#' @param theta a vector of values to be used as the copula dependence
#' parameter when computing the pdf. If theta is of length \eqn{< n} then the values
#' in the vector are repeated until the parameter vector is of length \eqn{n}
#'
#' @return The output is an \eqn{n x 1} vector of pdf values.
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
#' theta <- runif(n, min = 0, max = 10)
#' U <- matrix(runif(n * d), ncol = d)
#' pdfGumbel <- dGumbelCopula(U, theta)
dGumbelCopula <- function(U, theta) {
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

  if (any(theta < 1)) {
    stop('Parameter not in valid range')
  }

  if (length(theta) < n) {
    theta <- rep_len(theta, n)
  }

  pdfexpr <- GumbelPDFexpr(d)

  pdf_eval_fun <- function(i, U, theta, dim) {
    return(GumbelPDFeval(U[i, ], dim , theta[i]))
  }

  Y <- sapply(1:n, FUN = pdf_eval_fun, U = U, theta = theta, dim = d)
  return(Y)
}

#' GumbelLogLik
#'
#' A function to compute the Gumbel copula log-likelihood.
#'
#' This function computes the Gumbel copula log-likelihood at \eqn{n} points,
#' where the \eqn{i}th point in \eqn{d} dimensions is specified by the \eqn{i}th
#' row of the matrix \code{U}. The copula dependence parameter used to compute
#' the log-likelihood at each point is the \eqn{i}th value in the vector \code{theta}.
#' If \code{theta} is of length \eqn{< n}, values from \code{theta} are recycled
#' to produce a vector of length \eqn{n}
#'
#' @param theta a vector of copula parameter values each in the interval \eqn{[1, \infty)}.
#' If theta is of length \eqn{< n} then the values in theta will be recycled so
#' that the parameter vector has length \eqn{n}.
#' @param U an \eqn{n x d} matrix of points in \eqn{[0, 1]^d} at which to evaluate the
#' log-likelihood of the Gumbel copula.
#'
#' @return The output is an \eqn{n x 1} vector of values of the log-likelihood. Any
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
#' theta <- runif(n, min = 0, max = 10)
#' U <- matrix(runif(n * d), ncol = d)
#' loglikGumbel <- GumbelLogLik(theta, U)
GumbelLogLik <- function(theta, U) {
  # Check if U is a matrix and coerce to matrix form if necessary. We assume that
  # U is a (1 x d) vector in this instance.
  if (!any(class(U) %in% c('matrix'))) {
    U <- matrix(U, nrow = 1)
  }

  n <- nrow(U)
  d <- ncol(U)

  loglik <- log(dGumbelCopula(U, theta))

  loglik[is.infinite(loglik)] <- -1e9
  return(loglik)
}


# A helper function used to generate an expression for the cdf of the Gumbel copula
# for a given number of dimensions
GumbelCDFexpr <- function(dim) {
  CDFexpr <- paste0('(-log(u', 1:dim, '))^theta', collapse = ' + ')
  CDFexpr <- paste0('exp(-(', CDFexpr, ')^(1/theta))')
  CDFexpr <- parse(text = CDFexpr)
  return(CDFexpr)
}

# A helper function used to generate an expression for the pdf of the Gumbel copula
# for a given number of dimensions
GumbelPDFexpr <- function(dim) {
  CDFexpr <- GumbelCDFexpr(dim)
  PDFexpr <- CDFexpr
  for (i in 1:dim) {
    PDFexpr <- D(PDFexpr, paste0("u", i))
  }
  return(PDFexpr)
}

# A helper function used to evaluate the expression for the Gumbel copula pdf
# at a given point u with parameter theta
GumbelPDFeval <- function(u, dim, theta) {
  pdfexpr <- GumbelPDFexpr(dim)
  for (i in 1:dim) {
    assign(paste0("u", i), u[i])
  }
  y <- eval(pdfexpr)
  return(y)
}
