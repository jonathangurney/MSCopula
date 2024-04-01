RosenblattClayton <- function(U, theta) {
  n <- nrow(U)
  d <- ncol(U)

  if (d != 2) {
    stop('Rosenblatt transform only supported in bivariate case.')
  }

  if (length(theta) != n) {
    print('Dimensions of theta and U do not match.')
    theta <- rep_len(tbeta, length.out = n)
  }

  R <- matrix(nrow = n, ncol = d)
  R[ , 1] <- U[ , 1]

  R[ , 2] <- ((U[ , 1]^(-theta) + U[, 2]^(-theta) - 1)^(-(1/theta) - 1)) * U[ , 1]^(-theta - 1)
  return(R)
}
