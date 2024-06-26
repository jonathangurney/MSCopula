% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SJCCopulaDF.R
\name{pSJCCopula}
\alias{pSJCCopula}
\title{pSJCCopula}
\usage{
pSJCCopula(U, theta)
}
\arguments{
\item{U}{an \eqn{n x d} matrix giving n points in \eqn{[0, 1]^d} at which to evaluate
the cdf of the Symmetrized Joe-Clayton copula.}

\item{theta}{either a numeric vector of length 2 or an \eqn{n x 2} matrix of
parameter values in \eqn{[0, 1]}. Values in the first column of \code{theta}
determine the lower tail dependence, while those in the second column determine
the upper tail dependence. If \code{theta} is given as a matrix, values in the
\eqn{i}th row of \code{theta} are used to compute the cdf at the point specified
by row \eqn{i} of \code{U}. If \code{theta} is given as a numeric vector of
length 2 then these values are used to parameterise the cdf at all points.}
}
\value{
The output is an \eqn{n x 1} vector of cdf values.
}
\description{
A function to compute the distribution function of the Symmetrized Joe-Clayton copula.
}
\details{
This function computes the Symerized Joe-Clayton copula distribution function at \eqn{n}
points, where the \eqn{i}th point in \eqn{d} dimensions is specified by the
\eqn{i}th row of the matrix \code{U}. If \code{theta} is given as a numeric
vector then these parameter values are used to compute the cdf at all points.
If \code{theta} is a matrix, then the parameters in the \eqn{i}th row of
\code{theta} are used to compute the cdf at the point specified in the \eqn{i}th
row of \code{U}.
}
\examples{
set.seed(100)
n <- 10
d <- 3
theta <- matrix(runif(n * d), ncol = 2)
U <- matrix(runif(n * d), ncol = d)
cdfSymmetrizedJoeClayton <- pSJCCopula(U, theta)
}
\references{
\insertRef{nelsen2006}{MSCopula}
\insertRef{patton2006}{MSCopula}
}
