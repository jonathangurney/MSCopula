# Plotting and logging functions for MSCopula results

WriteLogs <- function(output, log_dir, data_start = as.Date(NA), data_end = as.Date(NA), indices = NA, freq = NA) {
  log_text <- character(length(output) - 4)
  switch(output$method,
         "EM" = {
           cats <- c('Method', 'Algorithm', 'Family', 'Start time', 'Runtime',
                     'Initial Params', 'Final Params', 'Log-likelihood', 'AIC', 'BIC', 'Regimes', 'Iterations',
                     'Convergence', 'Step tol', 'Data start', 'Data end', 'Indices', 'Frequency')
         },
         
         "Kim's Filter" = {
           cats <- c('Method', 'Algorithm', 'Family', 'Start time', 'Runtime',
                     'Initial Params', 'Final Params', 'Log-likelihood', 'AIC', 'BIC', 'Regimes', 'Function calls',
                     'Convergence', 'Grad tol', 'Step tol', 'Data start', 'Data end', 'Indices', 'Frequency')
         })
  
  for (i in 5:length(output)) {
    log_text[i - 4] <- paste0(output[[i]], collapse = ', ')
  }
  
  log_text <- append(log_text, format(data_start, '%Y-%m-%d'))
  log_text <- append(log_text, format(data_end, '%Y-%m-%d'))
  log_text <- append(log_text, paste0(indices, collapse = ', '))
  log_text <- append(log_text, freq)
  
  log_text <- paste0(cats, ': ', log_text, collapse = '\n')
  
  log_filename <- paste0(format(output[['init_time']], format = '%Y%m%d_%H%M%S'), '.txt')
  
  logConn <- file(paste0(log_dir, '/', log_filename))
  writeLines(log_text, logConn)
  close(logConn)
  
  # Determine number of parameters and regimes
  p <- length(output[['theta']])
  r <- output[['regimes']]
  
  theta_df <- output[['theta']][[1]]
  if (p > 1) {
    for (i in 2:p) {
      theta_df <- cbind(theta_df, output[['theta']][[i]])
    }
  }
  
  col_names <- character((3 + p) * r)
  for (i in 1:r) {
    col_names[i] <- paste0('eta', i)
    col_names[r + i] <- paste0('lambda', i)
    col_names[2*r + i] <- paste0('xi', i)
    for (j in 1:p) {
      col_names[3*r + (j - 1)*r + i] <- paste0('theta', j, '_', i)
    }
  }
  df <- as.data.frame(cbind(output[['eta']], output[['lambda']], output[['xi']], theta_df))
  colnames(df) <- col_names
  write.csv(df, file = paste0(log_dir, '/', format(output[['init_time']], format = '%Y%m%d_%H%M%S'), '.csv'), row.names = FALSE)
}

# ImportLogs <- function(log_dir) {
#   logs <- list.files(log_dir)
#   
#   for (filename in logs) {
#     filepath <- paste0(log_dir, '/', filename)
#     data <- read.table(filepath, sep = ': ')
#     print(data)
#   }
# }
