% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussianCopulaDF.R
\name{dGaussianCopula}
\alias{dGaussianCopula}
\title{dGaussianCopula}
\usage{
dGaussianCopula(U, rho)
}
\arguments{
\item{U}{an \eqn{n x d} matrix giving \eqn{n} points in \eqn{[0, 1]^d} at
which to evaluate the distribution function.}

\item{rho}{a \eqn{d x d} correlation matrix, a vector of length \eqn{n} or a
list of \eqn{n} \eqn{d x d} correlation matrices. \code{rho} can be a vector
of scalar values only if \eqn{d = 2}.}
}
\value{
The output is an \eqn{n x 1} vector of cdf values.
}
\description{
A function to compute the density function of the Gaussian copula.
}
\details{
This function computes the Gaussian copula density at \eqn{n} points, where
the \eqn{i}th point in \eqn{d} dimensions is specified by the \eqn{i}th row
of the matrix \code{U}. \code{rho} specifies the correlation matrices used to
parameterise the copula density function.
}
\examples{
set.seed(100)
n <- 10
d <- 3
rho <- replicate(n, p2rho(runif(d*(d - 1)/2, max = 0.6)), simplify = FALSE)
U <- matrix(runif(n * d), ncol = d)
pdfGaussian <- dGaussianCopula(U, rho)
}
\references{
\insertref{nelsen2006}{MSCopula}
}
