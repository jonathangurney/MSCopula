FitMSCopulaFilter <- function(U, pars_init, theta0, family, regimes = 2,
                              max_iter = 1e6, step_tol = 1e-10, grad_tol = 1e-10,
                              rel_tol = 1e-12, algorithm = 'Nelder-Mead',
                              log_output = TRUE) {
  # Initialize time
  start_time <- Sys.time()
  
  # Conduct minimization
  switch(algorithm,
         "nlm" = {
           optim_result <- stats::nlm(f = MSCopulaLogLik, p = pars_init, U = U,
                                      theta0 = theta0, family = family,
                                      regimes = regimes, gradtol = grad_tol,
                                      iterlim = max_iter, steptol = step_tol)
           
           pars_final <- optim_result$estimate
           fun_calls <- optim_result$iterations
           convergence <- optim_result$code
         },
         
         "Nelder-Mead" = {
           control <- list(maxit = max_iter, trace = 0, reltol = rel_tol)
           optim_result <- stats::optim(par = pars_init, fn = MSCopulaLogLik,
                                        method = 'Nelder-Mead', control = control,
                                        U = U, theta0 = theta0, family = family,
                                        regimes = regimes)
           
           pars_final <- optim_result$par
           fun_calls <- optim_result$counts[1]
           convergence <- optim_result$convergence
         },
         
         "BFGS" = {
           control <- list(maxit = max_iter, ndeps = rep(1e-7, length(pars_init)),
                           trace = 1)
           optim_result <- stats::optim(pars_init, MSCopulaLogLik, method = 'BFGS',
                                        control = control, U = U, theta0 = theta0,
                                        family = family, regimes = regimes)
           
           pars_final <- optim_result$par
           fun_calls <- optim_result$counts[1]
           convergence <- optim_result$convergence
         })
  
  fit_output <- MSCopulaLogLik(U, pars_final, theta0, family,
                               regimes = regimes, FALSE)
  
  runtime <- Sys.time() - start_time
  log_runtime <- difftime(Sys.time(), start_time, units = 'mins')
  
  # Calculate the AIC and BIC for the model
  aic <- 2 * length(pars_final) - 2 * fit_output$loglik
  bic <- log(nrow(U) * ncol(U)) * length(pars_final) - 2 * fit_output$loglik
  
  output <- list(eta = fit_output$eta,
               lambda = fit_output$lambda,
               xi = fit_output$xi,
               theta = fit_output$theta,
               method = "Kim's Filter",
               algorithm = algorithm,
               family = family,
               init_time = start_time,
               runtime = log_runtime,
               pars_init = pars_init,
               pars = pars_final,
               loglik = fit_output$loglik,
               aic = aic,
               bic = bic,
               regimes = regimes,
               fun_calls = fun_calls,
               convergence = convergence,
               grad_tol = grad_tol,
               step_tol = step_tol)
  
  return(output)
}