#' FitMSCopulaEM
#'
#' A function to fit a time-varying Markov-switching copula model using the
#' expectation maximisation algorithm.
#'
#' This function is used to fit a time-varying copula model where the dependence
#' parameter(s) for the copula vary according to a modified ARMA
#'
#' @param U (\eqn{n} x 2), a matrix of pseudo-observations in [0, 1]
#' @param pars_init vector of initial values for the parameters
#' @param theta0 vector of length \eqn{r} of initial values for the copula
#' dependence parameter in each regime
#' @param family the copula family to be used in the fitted model
#' @param regimes the number of regimes
#' @param max_iter the maximum number of iterations of the EM algorithm, the algorithm
#' will be terminated regardless of any conditions if this number is exceeded.
#' Default value is 1e6
#' @param max_fun_eval the maximum number of evaluations of the function \code{Q}.
#' Default value is 1e6
#' @param step_tol stopping criterion tolerance for convergence of the EM algorithm.
#' Default value is 1e-2
#' @param algorithm string determining the optimisation algorithm used to conduct
#' the minimisation in each step of the EM algorithm. Valid options are
#' \code{BFGS}, \code{Nelder-Mead} and \code{nlm}. Default is \code{Nelder-Mead}.
#' @param wall_time allocated wallclock time. Used if the algorithm is running on
#' a system where resource time is limited and checkpointing is required
#' @param save_time amount of time in hours to be set aside for saving data before
#' allocated wallclock time has expired
#'
#' @return \item{eta}{(\eqn{n} x \eqn{r}) matrix of filtered probabilities}
#' @return \item{lambda}{},
#' @return \item{xi = copula_loglik$xi,
#' @return \item{theta = EMresults$theta,
#' @return \item{method = 'EM',
#' @return \item{algorithm = algorithm,
#' @return \item{family = family,
#' @return \item{init_time = start_time,
#' @return \item{runtime = log_runtime,
#' @return \item{pars_init = pars_init,
#' @return \item{pars = pars_final,
#' @return \item{loglik = copula_loglik$loglik,
#' @return \item{aic = aic,
#' @return \item{bic = bic,
#' @return \item{regimes = regimes,
#' @return \item{iter = iter_counter,
#' @return \item{convergence = convergence,
#' @return \item{step_tol = step_tol
#'
#' @references
#' \insertRef{patton2006}{MSCopula}
#' \insertRef{silvafilho2012}{MSCopula}
#' \insertRef{nasri2020}{MSCopula}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @examples
FitMSCopulaEM <- function(U, pars_init, theta0, family, regimes = 2,
                          max_iter = 1e6, max_fun_eval = 1e6, step_tol = 1e-2,
                          algorithm = 'Nelder-Mead', wall_time = 1e10, save_time = 1) {
  # Initialize time to calculate runtime for algorithm
  start_time <- Sys.time()

  # Establish the dimensions of the problem
  n <- nrow(U)
  d <- ncol(U)
  r <- regimes

  # Track whether the algorithm converges or not
  convergence <- FALSE

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
             psi[i] <- mean(qnorm(U[(i - 10):(i - 1), 1]) * qnorm(U[(i - 10):(i - 1), 2]))
           }
         },

         "t" = {
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

  # Run the EM algorithm based on the optimization algorithm specified in the
  # function call

  # Initialize parameter values and function eval counter
  fun_eval_counter <- 0
  iter_counter <- 0
  pars <- pars_init

  for (i in 1:max_iter) {
    # Run the nlm algorithm

    print(paste0("iteration ", i))

    switch(algorithm,
           "nlm" = {
             optim_result <- stats::nlm(f = Q,
                                        p = pars,
                                        pars1 = pars,
                                        theta0 = theta0,
                                        U = U,
                                        psi = psi,
                                        family = family,
                                        regimes = regimes,
                                        optim_mode = TRUE)

             params <- optim_result$estimate
             min_value <- optim_result$minimum
             fun_eval_counter <- fun_eval_counter + optim_result$iterations
           },

           "Nelder-Mead" = {
             optim_result <- stats::optim(par = pars,
                                          fn = Q,
                                          method = 'Nelder-Mead',
                                          pars1 = pars,
                                          theta0 = theta0,
                                          U = U,
                                          psi = psi,
                                          family = family,
                                          regimes = regimes,
                                          optim_mode = TRUE)
             params <- optim_result$par
             min_value <- optim_result$value
             fun_eval_counter <- fun_eval_counter + optim_result$count[1]
           },

           "BFGS" = {
             optim_result <- stats::optim(par = pars,
                                          fn = Q,
                                          method = 'BFGS',
                                          pars1 = pars,
                                          theta0 = theta0,
                                          U = U,
                                          psi = psi,
                                          family = family,
                                          regimes = regimes,
                                          optim_mode = TRUE)
             params <- optim_result$par
             min_value <- optim_result$value
             fun_eval_counter <- fun_eval_counter + optim_result$count[1]
           })

    # Compute relevant values for convergence of EM algorithm
    step_norm <- sqrt(sum((pars - params)^2))

    iter_counter <- iter_counter + 1

    # Update parameter values
    pars <- params

    # Check if convergence has been achieved
    if (step_norm <= step_tol) {
      pars_final <- pars
      convergence <- TRUE
      break
    }

    # Check current runtime
    current_runtime <- difftime(Sys.time(), start_time, units = 'hours')

    # Conduct a smooth exit, saving necessary parmeters
    if (current_runtime > (wall_time - save_time)) {
      # We want to save the current parameter values and the runtime up to now
      pars_final <- pars
      break
    }
  }

  if (convergence) {
    print('EM Algorithm has converged')
  } else {
    message('Algorithm terminated to avoid overrunning wallclock time')
  }

  runtime <- Sys.time() - start_time
  log_runtime <- difftime(Sys.time(), start_time, units = 'mins')
  EMresults <- Q(pars_final, pars_final, U, theta0, regimes, psi, family, FALSE)

  copula_loglik <- MSCopulaLogLik(U, pars_final, theta0, regimes, family, FALSE)

  # Calculate the AIC and BIC for the model
  aic <- 2 * length(pars_final) - 2 * copula_loglik$loglik
  bic <- log(n * d) * length(pars_final) - 2 * copula_loglik$loglik

  output <- list(eta = EMresults$eta,
                 lambda = EMresults$lambda,
                 xi = copula_loglik$xi,
                 theta = EMresults$theta,
                 method = 'EM',
                 algorithm = algorithm,
                 family = family,
                 init_time = start_time,
                 runtime = log_runtime,
                 pars_init = pars_init,
                 pars = pars_final,
                 loglik = copula_loglik$loglik,
                 aic = aic,
                 bic = bic,
                 regimes = regimes,
                 iter = iter_counter,
                 convergence = convergence,
                 step_tol = step_tol)
  return(output)
}

# A function which computes the expected log-likelihood under the specification
# of par2 with respect to par1
Q <- function(pars1, pars2, U, theta0, regimes, psi, family, optim_mode) {
  n <- nrow(U)
  d <- ncol(U)
  r <- regimes

  # Specify number of time-varying parameters in the model
  if (family == 'sjc') {p <- 2}
  else {p <- 1}

  omega1 <- vector("list", p)
  alpha1 <- vector("list", p)
  beta1 <- vector("list", p)
  omega2 <- vector("list", p)
  alpha2 <- vector("list", p)
  beta2 <- vector("list", p)

  for (i in 1:p) {
    j <- (i - 1) * (r + 2) + 1
    omega1[[i]] <- pars1[j:(j + r - 1)]
    alpha1[[i]] <- pars1[(j + r)]
    beta1[[i]] <- pars1[(j + r + 1)]
    omega2[[i]] <- pars2[j:(j + r - 1)]
    alpha2[[i]] <- pars2[(j + r)]
    beta2[[i]] <- pars2[(j + r + 1)]
  }

  # To construct a valid transition matrix for r regimes, we need r * (r - 1)
  # coefficients
  P1 <- pars2matrix(pars1[((r + 2)*p + 1):((r + 2)*p + r * (r - 1))])
  P2 <- pars2matrix(pars2[((r + 2)*p + 1):((r + 2)*p + r * (r - 1))])

  # Create matrices of copula parameters, theta1 containing parameters for r
  # regimes with which to calculate the expectations and theta2 containing
  # parameters with which to compute the log-likelihood

  theta1 <- vector("list", p)
  theta2 <- vector("list", p)
  for (i in 1:p) {
    theta1[[i]] <- matrix(nrow = n, ncol = 2)
    theta2[[i]] <- matrix(nrow = n, ncol = 2)
    theta1[[i]][1, ] <- theta0
    theta2[[i]][1, ] <- theta0
  }

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
             return(1/(1.0001 + exp(-x)))
           }
         })

  for (j in 1:p) {
    for (i in 2:n) {
      lp1 <- omega1[[j]] + alpha1[[j]] * theta1[[j]][i - 1, ] + beta1[[j]] * psi[i]
      lp2 <- omega2[[j]] + alpha2[[j]] * theta2[[j]][i - 1, ] + beta2[[j]] * psi[i]
      theta1[[j]][i, ] <- transform_func(lp1)
      theta2[[j]][i, ] <- transform_func(lp2)
    }
  }

  # Matrix of copula densities
  c <- matrix(nrow = n, ncol = r)
  switch(family,
         "gaussian" = {
           for (i in 1:r) {
             c[ , i] <- dGaussianCopula(U, theta1[[1]][ , i])
           }
         },

         "gumbel" = {
           for (i in 1:r) {
             c[ , i] <- dGumbelCopula(U, theta1[[1]][ , i])
           }
         },

         "clayton" = {
           for (i in 1:r) {
             c[ , i] <- dClaytonCopula(U, theta1[[1]][ , i])
           }
         },

         "t" = {
           for (i in 1:r) {
             c[ , i] <- dtCopula(U, theta1[[1]][ , i], df = pars1[4*p + r * (r - 1) + i])
           }
         },

         "frank" = {
           for (i in 1:r) {
             c[ , i] <- dFrankCopula(U, theta1[[1]][ , i])
           }
         },

         "sjc" = {
           for (i in 1:r) {
             c[ , i] <- dSJCCopula(U, cbind(theta1[[1]][ , i], theta1[[2]][ , i]))
           }
         }
  )

  # Constructing the eta matrix
  eta <- matrix(nrow = n, ncol = r)
  eta0 <- rep(1/r, r)
  eta[1, ] <- c[1, ] * (t(P1) %*% eta0)/rep(sum(c[1, ] * (t(P1) %*% eta0)), 2)
  for (i in 2:n) {
    eta[i, ] <- (c[i, ] * (t(P1) %*% eta[i - 1, ]))/rep((t(c[i, ]) %*%
                                                    t(P1) %*% eta[i - 1, ]), 2)
  }

  # Constructing the eta_bar matrix
  eta_bar <- matrix(nrow = n, ncol = r)
  eta_bar[n, ] <- rep(1/r, r)
  for (i in 1:(n - 1)) {
    eta_bar[n - i, ] <- (c[n - i + 1, ] *
                           (P1 %*% eta_bar[n - i + 1, ]))/rep(sum(c[n - i + 1, ] *
                                                (P1 %*% eta_bar[n - i + 1, ])), r)
  }

  # Computing \lambda and \Lambda
  lambda <- (eta * eta_bar)/matrix(rep(rowSums(eta * eta_bar), 2), ncol = 2)
  Lambda <- vector("list", length = n)
  Lambda1_num <- P1 * (eta0 %*% t(eta_bar[1, ])) * matrix(rep(c[1, ], each = r),
                                                          ncol = r)
  Lambda[[1]] <- (Lambda1_num)/matrix(rep(sum(Lambda1_num), r^2), ncol = r)
  for (i in 2:n) {
    Lambda_num <- P1 * (eta[i - 1, ] %*% t(eta_bar[i, ])) * matrix(rep(c[i, ],
                                                            each = r), ncol = r)

    Lambda[[i]] <- (Lambda_num)/matrix(rep(sum(Lambda_num), r^2), ncol = r)
  }

  # Computing the log-likelihood
  LL <- matrix(nrow = n, ncol = r)
  switch(family,
         "gaussian" = {
           for (i in 1:r) {
             LL[ , i] <- GaussianLogLik(theta2[[1]][ , i], U)
           }
         },

         "gumbel" = {
           for (i in 1:r) {
             LL[ , i] <- GumbelLogLik(theta2[[1]][ , i], U)
           }
         },

         "clayton" = {
           for (i in 1:r) {
             LL[ , i] <- ClaytonLogLik(theta2[[1]][ , i], U)
           }
         },

         "t" = {
           for (i in 1:r) {
             LL[ , i] <- tLogLik(theta2[[1]][ , i], U,
                                 df = pars2[4*p + r * (r - 1) + i])
           }
         },

         "frank" = {
           for (i in 1:r) {
             LL[ , i] <- FrankLogLik(theta2[[1]][ , i], U)
           }
         },

         "sjc" = {
           for (i in 1:r) {
             LL[ , i] <- SJCLogLik(cbind(theta2[[1]][ , i],
                                         theta2[[2]][ , i]), U)
           }
         }
  )

  ExLL <- numeric(n)
  for (i in 1:n) {
    ExLL[i] <- sum(Lambda[[i]] * log(P2)) + sum(lambda[i, ] * LL[i, ])
  }

  if (optim_mode) {
    return(-sum(ExLL))
  }

  if (!optim_mode) {
    output <- list(eta = eta, lambda = lambda, theta = theta1, loglik = sum(ExLL))
    return(output)
  }
}
