RosenblattGaussian <- function(U, Rho) {
  n <- nrow(U)
  d <- ncol(U)

  if (d != 2) {
    stop('Rosenblatt transform only supported in bivariate case.')
  }

  if (length(Rho) != n) {
    print('Dimensions of rho and U do not match.')
    Rho <- rep_len(Rho, length.out = nrow(U))
  }

  R <- matrix(nrow = n, ncol = d)

  R[ , 1] <- U[ , 1]

  Q <- qnorm(U)
  R[ , 2] <- pnorm((Q[ , 2] - Q[ , 1] * Rho) * (1/sqrt(1 - Rho^2)))
  return(R)
}

