RosenblattGumbel <- function(U, theta) {
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

  V <- -log(U)
  W <- V[ , 1]^theta + V[ , 2]^theta

  R[ , 2] <- (1/U[ , 1]) * (W^(1/theta - 1)) * (V[ , 1]^(theta - 1)) * exp(-(W^(1/theta)))
  return(R)
}
