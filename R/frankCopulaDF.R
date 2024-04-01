#' pFrankCopula
#'
#' A function to compute the distribution function of the Frank copula.
#'
#' This function computes the Frank copula distribution function at \eqn{n}
#' points, where the \eqn{i}th point in \eqn{d} dimensions is specified by the
#' \eqn{i}th row of the matrix \code{U}. The copula dependence parameter used to
#' compute the distribution function at each point is the scalar \code{theta}.
#'
#' @param U an \eqn{n x d} matrix giving n points in \eqn{[0, 1]^d} at which to evaluate
#' the cdf of the Frank copula.
#' @param theta scalar dependence parameter for the Frank copula.
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
#' cdfFrank <- pFrankCopula(U, theta)
pFrankCopula <- function(U, theta) {
  # Check if U is a matrix and coerce to matrix form if necessary. We assume that
  # U is a (1 x d) vector in this instance.
  if (!any(class(U) %in% c("matrix"))) {
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

  X <- exp(-theta*U) - 1
  Y <- apply(X, FUN = prod, MARGIN = 1)
  Y <- - (1/theta) * log(1 + Y/((exp(-theta) - 1)^(d - 1)))
  return(Y)
}

#' dFrankCopula
#'
#' A function to compute the Frank copula density.
#'
#' This function computes the Frank copula density function at \eqn{n} points,
#' where the \eqn{i}th point in \eqn{d} dimensions is specified by the \eqn{i}th
#' row of the matrix \code{U}. The copula dependence parameter used to compute
#' the density at each point is the \eqn{i}th value in the vector \code{theta}.
#' If \code{theta} is of length \eqn{< n}, values from \code{theta} are recycled
#' to produce a vector of length \eqn{n}.
#'
#' @param U an \eqn{n x d} matrix giving n points in \eqn{[0, 1]^d} at which to evaluate
#' the pdf of the Frank copula
#' @param theta a vector of values to be used as the copula dependence
#' parameter when computing the pdf. If theta is of length \eqn{< n} then the values
#' in the vector are repeated until the parameter vector is of length \eqn{n}
#'
#' @return The output is an \eqn{n x 1} vector of pdf values. If any point returns
#' a value that evaluates to \code{Inf} then the function returns the value -1e9.
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
#' pdfFrank <- dFrankCopula(U, theta)
dFrankCopula <- function(U, theta) {
  n <- nrow(U)
  d <- ncol(U)

  if (d < 2) {
    stop('Copula must have at least 2 dimensions.')
  }

  if (length(theta) < n) {
    theta <- rep_len(theta, n)
  }

  pdf_eval_fun <- function(i, U, theta, dim) {
    return(FrankPDFeval(U[i, ], dim, theta[i]))
  }

  Y <- sapply(1:n, FUN = pdf_eval_fun, U = U, theta = theta, dim = d)
  return(Y)
}

#' FrankLogLik
#'
#' A function to compute the Frank copula log-likelihood.
#'
#' This function computes the Frank copula log-likelihood at \eqn{n} points,
#' where the \eqn{i}th point in \eqn{d} dimensions is specified by the \eqn{i}th
#' row of the matrix \code{U}. The copula dependence parameter used to compute
#' the log-likelihood at each point is the \eqn{i}th value in the vector \code{theta}.
#' If \code{theta} is of length \eqn{< n}, values from \code{theta} are recycled
#' to produce a vector of length \eqn{n}
#'
#' @param theta a vector of copula parameter values each in the interval \eqn{[0, \infty)}.
#' If theta is of length \eqn{< n} then the values in theta will be recycled so
#' that the parameter vector has length \eqn{n}.
#' @param U an \eqn{n x d} matrix of points in \eqn{[0, 1]^d} at which to evaluate the
#' log-likelihood of the Frank copula.
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
#' loglikFrank <- FrankLogLik(theta, U)
FrankLogLik <- function(theta, U) {
  n <- nrow(U)
  d <- ncol(U)

  loglik <- log(dFrankCopula(U, theta))

  loglik[is.infinite(loglik)] <- -1e9
  return(loglik)
}


FrankCDFexpr <- function(dim) {
  CDFexpr <- paste0('(exp(-theta * u', 1:dim, ') - 1)', collapse = ' * ')
  CDFexpr <- paste0('- (1/theta) * log(1 + (', CDFexpr, ')/((exp(-theta) - 1)^(', dim, ' - 1)))')
  CDFexpr <- parse(text = CDFexpr)
  return(CDFexpr)
}

FrankPDFexpr <- function(dim) {
  CDFexpr <- FrankCDFexpr(dim)
  PDFexpr <- CDFexpr
  for (i in 1:dim) {
    PDFexpr <- D(PDFexpr, paste0("u", i))
  }
  return(PDFexpr)
}

FrankPDFeval <- function(u, dim, theta) {
  pdfexpr <- FrankPDFexpr(dim)
  for (i in 1:dim) {
    assign(paste0("u", i), u[i])
  }
  y <- eval(pdfexpr)
  return(y)
}
