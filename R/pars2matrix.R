# Take vector of unconstrained parameters of length two and construct the
# corresponding transition matrix based on a logistic transformation
pars2matrix <- function(x) {
  r <- (1 + sqrt(4*length(x) + 1))/2
  logistic <- function(x) {
    return(1/(1.00001 + exp(-x)))
  }
  p <- logistic(x)
  P <- matrix(p, ncol = r - 1, byrow = T)
  P <- cbind(P, 1 - rowSums(P))
  if (any((P > 1) | (P < 0))) {
    stop('Invalid transition matrix')
  }
  return(P)
}
