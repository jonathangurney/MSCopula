MSCopulaLogLik <- function(U, pars, theta0, regimes, family, optim_mode = TRUE) {
  n <- nrow(U)
  d <- ncol(U)
  r <- regimes
  
  #Specify the number of time-varying parameters
  if (family %in% c('sjc')) {p <- 2}
  else {p <- 1}
  
  omega <- vector("list", p)
  alpha <- vector("list", p)
  beta <- vector("list", p)
  
  for (i in 1:p) {
    j <- (i - 1)*(r + 2) + 1
    omega[[i]] <- pars[j:(j + r - 1)]
    alpha[[i]] <- pars[(j + r)]
    beta[[i]] <- pars[(j + r + 1)] 
  }
  
  P <- pars2matrix(pars[(p*(r + 2) + 1):(p*(r + 2) + r * (r - 1))])
  
  switch(family,
         "t" = {
           df <- pars[(p*(r + 2) + r * (r - 1) + 1):(p*(r + 2) + r^2)]
         })
  
  # Compute the forcing term vector psi. The method to compute this changes
  # based on the copula family, see Patton (2006).
  
  archim_cops <- c('clayton', 'frank', 'gumbel', 'sjc')
  
  if (family %in% archim_cops) {
    psi <- numeric(n)
    psi[1] <- 0
    for (i in 2:n) {
      if (i <= 10) {
        psi[i] <- mean(abs(U[1:(i - 1), 1] - U[1:(i - 1), 2]))
        next
      }
      psi[i] <- mean(abs(U[(i - 10):(i - 1), 1] - U[(i - 10):(i - 1), 2]))
    }
  }
  
  switch(family,
         "gaussian" = {
           # The transformation function to ensure the parameter is constrained
           # to a viable range for the copula family.
           psi <- numeric(n)
           psi[1] <- 0
           for (i in 2:n) {
             if (i <= 10) {
               psi[i] <- mean(qnorm(U[1:(i - 1), 1]) * qnorm(U[1:(i - 1), 2]))
               next
             }
             psi[i] <- mean(qnorm(U[(i - 10):(i - 1), 1]) *
                              qnorm(U[(i - 10):(i - 1), 2]))
           }
         },
         
         "t" = {
           # The use of df = 4 in the following is really a placeholder and
           # in the literature (de Silva Filho et at, 2013) is replaced by the
           # number of degrees of freedom from the constant t-copula fitted to
           # the data.
           
           psi <- numeric(n)
           psi[1] <- 0
           for (i in 2:n) {
             if (i <= 10) {
               psi[i] <- mean(qt(U[1:(i - 1), 1], df = 4) *
                                qt(U[1:(i - 1), 2], df = 4))
               next
             }
             psi[i] <- mean(qt(U[(i - 10):(i - 1), 1], df = 4) *
                              qt(U[(i - 10):(i - 1), 2], df = 4))
           }
         }
  )
  
  # Specify the transformation used to constrain the parameter.
  
  switch(family, 
         "gaussian" = {
           transform_func <- function(x) {
             return(1.9998/(exp(-x) + 1) - 0.9999)
           }
         },
         
         "gumbel" = {
           transform_func <- function(x) {
             return(1 + x^2)
           }
         },
         
         "clayton" = {
           transform_func <- function(x) {
             return(x^2)
           }
         },
         
         "t" = {
           transform_func <- function(x) {
             return(1.9998/(exp(-x) + 1) - 0.9999)
           }
         },
         
         "frank" = {
           transform_func <- function(x) {
             return(x)
           }
         },
         
         "sjc" = {
           transform_func <- function(x) {
             return(0.9998/(1 + exp(-x)) + 0.0001)
           }
         })
  
  theta <- vector("list", p)
  for (i in 1:p) {
    theta[[i]] <- matrix(nrow = n, ncol = 2)
    theta[[i]][1, ] <- theta0
    for (j in 2:n) {
      arma <- omega[[i]] + alpha[[i]] * theta[[i]][(j - 1), ] + beta[[i]] * psi[j]
      theta[[i]][j, ] <- transform_func(arma)
    }
  }
  
  # Matrix of copula densities
  c <- matrix(nrow = n, ncol = r)
  switch(family,
         "gaussian" = {
           for (i in 1:r) {
             c[ , i] <- dGaussianCopula(U, theta[[1]][ , i])
           }
         },
         
         "gumbel" = {
           for (i in 1:r) {
             c[ , i] <- dGumbelCopula(U, theta[[1]][ , i])
           }
         },
         
         "clayton" = {
           for (i in 1:r) {
             c[ , i] <- dClaytonCopula(U, theta[[1]][ , i])
           }
         },
         
         "t" = {
           for (i in 1:r) {
             c[ , i] <- dtCopula(U, theta[[1]][ , i],
                                 df = pars[(r + 2)*p + r * (r - 1) + i])
           }
         },
         
         "frank" = {
           for (i in 1:r) {
             c[ , i] <- dFrankCopula(U, theta[[1]][ , i])
           }
         },
         
         "sjc" = {
           for (i in 1:r) {
             c[ , i] <- dSJCCopula(U, cbind(theta[[1]][ , i], theta[[2]][ , i]))
           }
         }
  )
  
  # Constructing the eta matrix and matrix of predictive probabilities
  eta <- matrix(nrow = n, ncol = r)
  xi <- matrix(nrow = n, ncol = r)
  
  eta0 <- rep(1/r, r)
  
  eta[1, ] <- c[1, ] * (t(P) %*% eta0)/rep(sum(c[1, ] * (t(P) %*% eta0)), r)
  xi[1, ] <- t(P) %*% eta0
  for (i in 2:n) {
    eta[i, ] <- (c[i, ] * (t(P) %*% eta[i - 1, ]))/rep((t(c[i, ]) %*%
                                                      t(P) %*% eta[i - 1, ]), r)
    xi[i, ] <- t(P) %*% eta[i - 1, ]
  }
  
  # Constructing the eta_bar matrix
  eta_bar <- matrix(nrow = n, ncol = r)
  eta_bar[n, ] <- rep(1/r, r)
  for (i in 1:(n - 1)) {
    eta_bar[n - i, ] <- (c[n - i + 1, ] *
                           (P %*% eta_bar[n - i + 1, ]))/rep(sum(c[n - i + 1, ] *
                           (P %*% eta_bar[n - i + 1, ])), r)
  }
  
  # Constructing the matrix of smoothed probabilities
  lambda <- matrix(nrow = n, ncol = r)
  lambda[n, ] <- eta[n, ]
  for (i in 2:n) {
    lambda[n - i + 1, ] <- eta[n - i + 1, ] *
      (P %*% (lambda[n - i + 2, ]/(t(P) %*% eta[n - i + 1, ])))
  }
  
  # Constructing the log-likelihood
  LogLik <- log(rowSums(c * xi))
  
  if (optim_mode) {
    # Return the negative log-likelihood for minimization
    return(-sum(LogLik))
  }
  
  if (!optim_mode) {
    # Return relevant probabilities and parameter dynamics
    output <- list(eta = eta, xi = xi, lambda = lambda,
                   theta = theta, pars = pars, loglik = sum(LogLik))
  }
}

