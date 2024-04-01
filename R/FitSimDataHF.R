# Define function that will:
#   1. Simulate a data set using given parameters, a given number of observations
#      and a given copula family
#   2. Fit a MS copula model via the Hamilton filter
#   3. Fit a MS copula model via the EM algorithm
#   4. Return the following:
#       - Mean dependence parameter error for both models
#       - Mean state error for each model
#       - Final parameter values for each model

FitSimData <- function(n, pars, pars_init, trans_mat, theta0, family, regimes) {
  # Simulate the data
  data <- SimData(n = n, pars = pars, trans_mat = trans_mat, family = family)
  
  # Fit the MS copula model using HF
  HFResult <- FitMSCopulaFilter(U = data$U,
                                pars_init = pars_init,
                                theta0 = theta0,
                                family = family,
                                regimes = regimes)
  
  # Processing the data for output
  
  # Find the state with highest probability at each time period
  HF_states <- apply(X = HFResult$eta, MARGIN = 1, FUN = function(x) {sample(which(x == max(x)), 1)})
  
  # Find the dependence parameter corresponding to the state with the highest
  # probability at each time  
  HF_theta <- sapply(X = 1:n, FUN = function(x) {return(HFResult$theta[[1]][x, HF_states[x]])})
  
  # Find the mean error, RMSE, MAE
  HF_ME <- mean(HF_theta - data$theta_sim)
  HF_RMSE <- sqrt(mean((HF_theta - data$theta_sim)^2))
  HF_MAE <- mean(abs(HF_theta - data$theta_sim))
  
  # Create list for output
  output <- list(HF_errors = c(HF_ME, HF_RMSE, HF_MAE),
                 HF_pars = HFResult$pars,
                 runtimes = HFResult$runtime)
  
  return(output)
}


