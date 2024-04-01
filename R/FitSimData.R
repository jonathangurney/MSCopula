# A function that will:
#   1. Simulate a data set using given parameters, a given number of observations
#      and a given copula family
#   2. Fit a MS copula model via the Hamilton filter
#   3. Fit a MS copula model via the EM algorithm
#   4. Return the following:
#       - Mean dependence parameter error for both models
#       - Mean state error for each model
#       - Final parameter values for each model

FitSimData <- function(n, pars, pars_init, trans_mat, theta0, family, regimes, ...) {
  # Set default options for keyword parameters
  
  # Declare appropriate keyword arguments
  key_args <- c("HF_optim_algo", "HF_NM_rel_tol", "EM_optim_algo",
                "EM_NM_rel_tol", "EM_step_tol")
  
  # Set default values for keyword arguments
  HF_optim_algo <- "Nelder-Mead"
  HF_NM_rel_tol <- 1e-15
  EM_optim_algo <- "Nelder-Mead"
  EM_NM_rel_tol <- 1e-9
  EM_step_tol <- 1e-3
  
  key_params <- (...)
  for (name in names(key_params)) {
    if (!(name %in% key_args)) {
      stop(paste0('Argument ', name, ' not accepted by FitSimData'))
    }
    assign(name, key_params[[name]])
  }
  
  data <- SimData(n = n, pars = pars, trans_mat = trans_mat,
                  family = family, regimes = regimes)
  
  # Fit the MS Copula model using the HF algorithm
  HFResult <- FitMSCopulaFilter(U = data$U,
                                pars_init = pars_init,
                                theta0 = theta0,
                                family = family,
                                regimes = regimes,
                                rel_tol = HF_NM_rel_tol)
  
  # Fit the MS copula model using the EM algorithm
  EMResult <- FitMSCopulaEM(U = data$U,
                            pars_init = pars_init,
                            theta0 = theta0,
                            family = family,
                            regimes = regimes,
                            max_iter = 1e5,
                            algorithm = EM_optim_algo,
                            step_tol = EM_step_tol)

  # Processing the data for output

  # Find the state with highest probability at each time period
  HF_states <- apply(X = HFResult$eta, MARGIN = 1,
                     FUN = function(x) {sample(which(x == max(x)), 1)})
  EM_states <- apply(X = EMResult$eta, MARGIN = 1,
                     FUN = function(x) {sample(which(x == max(x)), 1)})

  # Find the dependence parameter corresponding to the state with the highest
  # probability at each time
  HF_theta <- sapply(X = 1:n,
                     FUN = function(x) {return(HFResult$theta[[1]][x, HF_states[x]])})
  EM_theta <- sapply(X = 1:n,
                     FUN = function(x) {return(EMResult$theta[[1]][x, EM_states[x]])})

  # Find the mean error, RMSE, MAE
  HF_ME <- mean(HF_theta - data$theta_sim)
  EM_ME <- mean(EM_theta - data$theta_sim)
  HF_RMSE <- sqrt(mean((HF_theta - data$theta_sim)^2))
  EM_RMSE <- sqrt(mean((EM_theta - data$theta_sim)^2))
  HF_MAE <- mean(abs(HF_theta - data$theta_sim))
  EM_MAE <- mean(abs(EM_theta - data$theta_sim))
  
  # Find the log-likelihood, AIC, and BIC for the fitted models
  HF_loglik <- HFResult$loglik
  EM_loglik <- EMResult$loglik
  HF_AIC <- HFResult$aic
  EM_AIC <- EMResult$aic
  HF_BIC <- HFResult$bic
  EM_BIC <- EMResult$bic

  # Create list for output
  output <- list(HF_errors = c(HF_ME, HF_RMSE, HF_MAE),
                 EM_errors = c(EM_ME, EM_RMSE, EM_MAE),
                 HF_pars = HFResult$pars,
                 EM_pars = EMResult$pars,
                 runtimes = c(HFResult$runtime, EMResult$runtime),
                 HF_GOF = c(HF_loglik, HF_AIC, HF_BIC),
                 EM_GOF = c(EM_loglik, EM_AIC, EM_BIC))
  
  return(output)
}


