# Create figures from EM results

WriteFigs <- function(output, fig_dir, index1, index2) {
  r <- output[['regimes']]
  p <- length(output[['theta']])
  
  # Set graphical parameters
  linewidth <- '1.5'
  
  # Create plot of filtered probabilities
  fp_filename <- paste0(output[['family']], '_Filtered_Probabilities.png')
  png(filename = paste0(fig_dir, '/', fp_filename), width= 1200, height = 800, pointsize = 16)
  par(mfrow = c(r, 1))
  for (i in 1:r) {
    main <- paste0('State ', i - 1, ' - Filtered probabilities for ', index1, ' and ', index2, ' with ', output[['family']], ' copula')
    plot(dates, output[['eta']][ , i], type = 'l', ylim = c(0, 1), ylab = 'Probability', lwd = linewidth,
         main = main)
  }
  dev.off()
  
  # Create plot of smoothed probabilities
  sp_filename <- paste0(output[['family']], '_Smoothed_Probabilities.png')
  png(filename = paste0(fig_dir, '/', sp_filename), width = 1200, height = 800, pointsize = 16)
  par(mfrow = c(r, 1))
  for (i in 1:r) {
    main <- paste0('State ', i - 1, ' - Smoothed probabilities for ', index1, ' and ', index2, ' with ', output[['family']], ' copula')
    plot(dates, output[['lambda']][ , i], type = 'l', ylim = c(0, 1), ylab = 'Probability', lwd = linewidth,
         main = main)
  }
  dev.off()
  
  # Set y limits for the dependence parameters
  switch(output[['family']],
         "gaussian" = {
           par_ylim <- vector("list", r)
           for (i in 1:r) {
             par_ylim[[i]] <- c(-1, 1)
           }
         },
         
         "gumbel" = {
           par_ylim <- vector("list", r)
           for (i in 1:r) {
             theta_min <- min(output[['theta']][[1]][ , i])
             theta_max <- max(output[['theta']][[1]][ , i])
             theta_diff <- theta_max - theta_min
             par_ylim[[i]] <- c(theta_min - 0.1 * theta_diff,  theta_max + 0.1 * theta_diff)
           }
         },
         
         "clayton" = {
           par_ylim <- vector("list", r)
           for (i in 1:r) {
             theta_min <- min(output[['theta']][[1]][ , i])
             theta_max <- max(output[['theta']][[1]][ , i])
             theta_diff <- theta_max - theta_min
             par_ylim[[i]] <- c(theta_min - 0.1 * theta_diff,  theta_max + 0.1 * theta_diff)
           }
         },
         
         "frank" = {
           par_ylim <- vector("list", r)
           for (i in 1:r) {
             theta_min <- min(output[['theta']][[1]][ , i])
             theta_max <- max(output[['theta']][[1]][ , i])
             theta_diff <- theta_max - theta_min
             par_ylim[[i]] <- c(theta_min - 0.1 * theta_diff,  theta_max + 0.1 * theta_diff)
           }
         },
         
         "t" = {
           par_ylim <- vector("list", r)
           for (i in 1:r) {
             par_ylim[[i]] <- c(-1, 1)
           }
         },
         
         "sjc" = {
           par_ylim <- vector("list", r)
           for (i in 1:r) {
             par_ylim[[i]] <- c(0, 1)
           }
         })
  
  for (i in 1:p) {
    if (p != 1) {
      par_filename <- paste0(output[['family']], '_Parameter', i, '.png')
    }
    else {
      par_filename <- paste0(output[['family']], '_Parameter.png')
    }
    
    png(filename = paste0(fig_dir, '/', par_filename), width= 1200, height = 800, pointsize = 16)
    par(mfrow = c(r, 1))
    
    for (j in 1:r) {
      if (p != 1) {
        main <- paste0('State ', j - 1, ' - Dependence parameter ', i, ' for ', index1, ' and ', index2, ' with ', output[['family']], ' copula')
      }
      else {
        main <- paste0('State ', j - 1, ' - Dependence parameter for ', index1, ' and ', index2, ' with ', output[['family']], ' copula')
      }
      plot(dates, output[['theta']][[i]][ , j], type = 'l', ylim = par_ylim[[i]], xlab = 'Date', ylab = 'Parameter value', lwd = linewidth,
           main = main)
    }
    dev.off()
  }
}