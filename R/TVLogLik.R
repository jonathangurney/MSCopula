EstTVCopula <- function(U, params0, theta0, family) {
  n <- nrow(U)
  d <- ncol(U)

  # Create vector for the forcing terms
  switch(family,
         "gaussian" = {
           psi <- numeric(n)
           psi[1] <- 0
           for (i in 2:n) {
             if (i <= 10) {
               psi[i] <- mean(qnorm(U[1:(i - 1), 1]) * qnorm(U[1:(i - 1), 2]))
               next
             }
             psi[i] <- mean(qnorm(U[(i - 10):(i - 1), 1]) * qnorm(U[(i - 10):(i - 1), 2]))
           }

           transform.func <- function(x) {
             return(2/(exp(-x) + 1) - 1)
           }
         },

         "gumbel" = {
           psi <- numeric(n)
           for (i in 1:nrow(U)) {
             if (i == 1) {
               psi[i] <- 0
               next
             }

             if (i <= 10) {
               psi[i] <- mean(abs(U[1:(i - 1), 1] - U[1:(i - 1), 2]))
               next
             }

             psi[i] <- mean(abs(U[(i - 10):(i - 1), 1] - U[(i - 10):(i - 1), 2]))
           }

           transform.func <- function(x) {
             return(1 + x^2)
           }
         })

  # Construct function that evaluates the log-likelihood

  LLeval <- function(params) {
    omega <- params[1]
    alpha <- params[2]
    beta <- params[3]

    theta <- numeric(n)

    theta[1] <- theta0
    for (i in 2:n) {
      lp <- omega + alpha * theta[i - 1] + beta * psi[i]
      theta[i] <- transform.func(lp)
    }

    switch(family,
           "gaussian" = {
             LL <- GaussianLogLik(theta, U)
           },

           "gumbel" = {
             LL <- GumbelLogLik(theta, U)
           })
    # print(paste0('w = ', params[1], ', a = ', params[2], ', b = ', params[3]))
    return(-sum(LL))
  }

  # Now use stats::nlm or stats::optim to optimize the function LLeval

  lower <- rep(c(-25), 3)
  upper <- rep(c(25), 3)

  result <- stats::nlm(f = LLeval, p = params0)
  return(result)
}
