library(quantmod)
library(stats)

start_date <- as.Date('01-01-2010', format = "%d-%m-%Y")
end_date <- as.Date('31-12-2019', format = "%d-%m-%Y")

getSymbols(Symbols = c('^FTSE', '^BVSP', '^GSPC'),
           src = 'yahoo',
           from = start_date,
           to = end_date,
           periodicity = 'daily')

# Check if data contains any NA values

if (any(is.na(GSPC))) {
  GSPC <- GSPC[complete.cases(GSPC), ]
}

if (any(is.na(FTSE))) {
  FTSE <- FTSE[complete.cases(FTSE), ]
}

if (length(setdiff(zoo::index(FTSE), zoo::index(GSPC))) != 0) {
  common_dates <- as.Date(intersect(zoo::index(FTSE), zoo::index(GSPC)),
                          origin = "1970-01-01")
  GSPC <- GSPC[common_dates, ]
  FTSE <- FTSE[common_dates, ]
}

# Create vector of log close prices
GSPC$GSPC.LogClose <- log(GSPC$GSPC.Close)
FTSE$FTSE.LogClose <- log(FTSE$FTSE.Close)

# Create vector of log returns
GSPC$GSPC.LogReturn <- c(diff(GSPC$GSPC.LogClose))
FTSE$FTSE.LogReturn <- c(diff(FTSE$FTSE.LogClose))

GSPCLogRet <- GSPC$GSPC.LogReturn[2:length(FTSE$FTSE.LogReturn)]
FTSELogRet <- FTSE$FTSE.LogReturn[2:length(FTSE$FTSE.LogReturn)]

n <- length(FTSELogRet)

# Plot the acf and pacf for the returns and the squared returns
max_lag <- 30
par(mfrow = c(1, 2))
acf(GSPCLogRet, col = 'green', lwd = 2, lag.max = max_lag)
pacf(GSPCLogRet, col = 'red', lwd = 2, lag.max = max_lag)
acf(FTSELogRet, col = 'green', lwd = 2, lag.max = max_lag)
pacf(FTSELogRet, col = 'red', lwd = 2, lag.max = max_lag)
acf(GSPCLogRet^2, col = 'green', lwd = 2, lag.max = max_lag)
pacf(GSPCLogRet^2, col = 'red', lwd = 2, lag.max = max_lag)
acf(FTSELogRet^2, col = 'green', lwd = 2, lag.max = max_lag)
pacf(FTSELogRet^2, col = 'red', lwd = 2, lag.max = max_lag)

# Models to consider for the FTSE
FTSE_ARMA_orders <- list(ar = c(0), ma = c(0))
FTSE_GARCH_orders <- list(ar = c(1), ma = c(1))

# Get indices
n_FTSE_p <- length(FTSE_ARMA_orders$ar)
n_FTSE_q <- length(FTSE_ARMA_orders$ma)
n_FTSE_r <- length(FTSE_GARCH_orders$ar)
n_FTSE_s <- length(FTSE_GARCH_orders$ma)

# Specify conditional distributions to be used
marginal_dists <- c('norm', 'snorm', 'ged', 'sged', 'sstd', 'std')
n_dists <- length(marginal_dists)

# Create vectors useful for creating array dimensions

FTSE_array_dims <- c(length(FTSE_ARMA_orders$ar),
                     length(FTSE_ARMA_orders$ma),
                     length(FTSE_GARCH_orders$ar),
                     length(FTSE_GARCH_orders$ma),
                     length(marginal_dists))

# Initialize array for storing residuals
FTSE_residuals <- array(dim = c(FTSE_array_dims, n))

# Initialize array for storing information criteria statistics
FTSE_IC <- array(dim = c(FTSE_array_dims, 2))

# Initialize array for storing Ljung-Box p-values
FTSE_LB <- array(dim = c(FTSE_array_dims))

# Initialize array for storing Kolmogorov-Smirnov and Anderson-Darling p-values
FTSE_KS_AD <- array(dim = c(FTSE_array_dims, 2))

# Initialize array for storing ARCH-LM p-values
FTSE_ARCH_LM <- array(dim = FTSE_array_dims)

# Initialize array for storing formulae
FTSE_formulae <- array(dim = FTSE_array_dims)

# Now look to fit ARMA-GARCH models to the data, we fit ARMA(p, q) and GARCH(r, s)
for (i in 1:n_dists) {
  print(marginal_dists[i])
  for (j in 1:n_FTSE_p) {
    for (k in 1:n_FTSE_q) {
      for (l in 1:n_FTSE_r) {
        for (m in 1:n_FTSE_s) {
          err <- FALSE
          
          # Set the orders for the current model
          p <- FTSE_ARMA_orders$ar[j]
          q <- FTSE_ARMA_orders$ma[k]
          r <- FTSE_GARCH_orders$ar[l]
          s <- FTSE_GARCH_orders$ma[m]
          
          # Create the formula for the model to be fit
          formula_text <- paste0('~ arma(', p, ', ', q, ') + garch(', r, ', ', s, ')')
          form <- formula(formula_text)
          FTSE_formulae[j, k, l, m, i] <- formula_text
          print(formula_text)
          
          # Fit the model described above using garchFit
          tryCatch(
            FTSE_model <- fGarch::garchFit(formula = form,
                                           data = FTSELogRet,
                                           trace = F,
                                           cond.dist = marginal_dists[i]),
            
            error = function(e) {
              err <- TRUE
              message('An error occurred:\n', e)
            }
          )
          
          if (!err) {
            #
            # Store output and statistics arising from model fitting
            #
            
            # Store vector of residuals
            FTSE_res <- fGarch::residuals(FTSE_model, standardize = TRUE)
            FTSE_residuals[j, k, l, m, i, ] <- FTSE_res
            
            # Store information criteria
            FTSE_IC[j, k, l, m, i, ] <- c(FTSE_model@fit$ics[['AIC']],
                                          FTSE_model@fit$ics[['BIC']])
            
            # Compute Ljung-Box p-value for residuals
            FTSE_LBtest10 <- stats::Box.test(FTSE_res, lag = 10,
                                             type = 'Ljung-Box', fitdf = (p + q))
            
            #Record p-values from the Ljung-Box test
            FTSE_LB[j, k, l, m, i] <- FTSE_LBtest10$p.value
            
            # Compute ARCH-LM statistic for the residuals
            FTSE_ARCH_LM[j, k, l, m, i] <- fDMA::archtest(FTSE_res, lag = 1)$p.value
            
            # Compute the Kolmogorov-Smirnov and Anderson-Darling statistics
            
            # Use the probability integral transform to transform to uniform random
            # variables
            
            # Use PIT to create uniform variables
            
            # For the 'norm' conditional distribution
            if (marginal_dists[i] == 'norm') {
              U <- pnorm(FTSE_res, mean=0, sd=1)
            }
            
            # For the 'snorm' conditional distribution
            if (marginal_dists[i] == 'snorm') {
              skew <- FTSE_model@fit$coef[['skew']]
              U <- fGarch::psnorm(FTSE_res, mean=0, sd=1, xi=skew)
            }
            
            # For the 'ged' conditional distribution
            if (marginal_dists[i] == 'ged') {
              shape <- FTSE_model@fit$coef[['shape']]
              U <- fGarch::pged(FTSE_res, mean=0, sd=1, nu=shape)
            }
            
            # For the 'sged' conditional distribution
            if (marginal_dists[i] == 'sged') {
              shape <- FTSE_model@fit$coef[['shape']]
              skew <- FTSE_model@fit$coef[['skew']]
              U <- fGarch::psged(FTSE_res, mean=0, sd=1, nu=shape, xi=skew)
            }
            
            # For the 'sstd' conditional distribution
            if (marginal_dists[i] == 'sstd') {
              shape <- FTSE_model@fit$coef[['shape']]
              skew <- FTSE_model@fit$coef[['skew']]
              U <- fGarch::psstd(FTSE_res, mean=0, sd=1, nu=shape, xi=skew)
            }
            
            # For the 'std' conditional distribution
            if (marginal_dists[i] == 'std') {
              shape <- FTSE_model@fit$coef[['shape']]
              U <- fGarch::pstd(FTSE_res, mean=0, sd=1, nu=shape)
            }
            
            # Conduct the KS test and record p-values
            FTSE_KS_AD[j, k, l, m, i, 1] <- KScorrect::LcKS(U, cdf='punif')$p.value
            
            # Conduct the AD test and record p-values
            FTSE_KS_AD[j, k, l, m, i, 2] <- ADGofTest::ad.test(U, null='punif')$p.value
          }
        }
      }
    }
  }
}


#############################################################################
#
# Same but now it's GSPC
#
#############################################################################


GSPC_ARMA_orders <- list(ar = c(0, 1), ma = c(0))
GSPC_GARCH_orders <- list(ar = c(1), ma = c(1))

# Get indices
n_GSPC_p <- length(GSPC_ARMA_orders$ar)
n_GSPC_q <- length(GSPC_ARMA_orders$ma)
n_GSPC_r <- length(GSPC_GARCH_orders$ar)
n_GSPC_s <- length(GSPC_GARCH_orders$ma)

# Initialize arrays for storing GOF test results

marginal_dists <- c('norm', 'snorm', 'ged', 'sged', 'sstd', 'std')
n_dists <- length(marginal_dists)

# Create vectors useful for creating array dimensions

GSPC_array_dims <- c(length(GSPC_ARMA_orders$ar),
                     length(GSPC_ARMA_orders$ma),
                     length(GSPC_GARCH_orders$ar),
                     length(GSPC_GARCH_orders$ma),
                     length(marginal_dists))

# Initialize array for storing residuals
GSPC_residuals <- array(dim = c(GSPC_array_dims, n))

# Initialize array for storing information criteria statistics
GSPC_IC <- array(dim = c(GSPC_array_dims, 2))

# Initialize array for storing Ljung-Box p-values
GSPC_LB <- array(dim = c(GSPC_array_dims))

# Initialize array for storing ARCH-LM p-values
GSPC_ARCH_LM <- array(dim = GSPC_array_dims)

# Initialize array for storing Kolmogorov-Smirnov and Anderson-Darling p-values
GSPC_KS_AD <- array(dim = c(GSPC_array_dims, 2))

# Initialize array for storing formulae
GSPC_formulae <- array(dim = GSPC_array_dims)

# Now look to fit ARMA-GARCH models to the data, we fit ARMA(p, q) and GARCH(r, s)
for (i in 1:n_dists) {
  print(marginal_dists[i])
  for (j in 1:n_GSPC_p) {
    for (k in 1:n_GSPC_q) {
      for (l in 1:n_GSPC_r) {
        for (m in 1:n_GSPC_s) {
          err <- FALSE
          
          # Set the orders for the current model
          p <- GSPC_ARMA_orders$ar[j]
          q <- GSPC_ARMA_orders$ma[k]
          r <- GSPC_GARCH_orders$ar[l]
          s <- GSPC_GARCH_orders$ma[m]
          
          # Create the formula for the model to be fit
          formula_text <- paste0('~ arma(', p, ', ', q, ') + garch(', r, ', ', s, ')')
          form <- formula(formula_text)
          GSPC_formulae[j, k, l, m, i] <- formula_text
          print(formula_text)
          
          # Fit the model described above using garchFit
          tryCatch(
            GSPC_model <- fGarch::garchFit(formula = form,
                                           data = GSPCLogRet,
                                           trace = F,
                                           cond.dist = marginal_dists[i]),
            
            error = function(e) {
              err <- TRUE
              message('An error occurred:\n', e)
            }
          )
          
          if (!err) {
            #
            # Store output and statistics arising from model fitting
            #
            
            # Store vector of residuals
            GSPC_res <- fGarch::residuals(GSPC_model, standardize = TRUE)
            GSPC_residuals[j, k, l, m, i, ] <- GSPC_res
            
            # Store information criteria
            GSPC_IC[j, k, l, m, i, ] <- c(GSPC_model@fit$ics[['AIC']],
                                          GSPC_model@fit$ics[['BIC']])
            
            # Compute Ljung-Box p-value for residuals and squared residuals
            GSPC_LBtest10 <- stats::Box.test(GSPC_res, lag = 10,
                                             type = 'Ljung-Box', fitdf = (p + q))
            
            #Record p-values fromt the Ljung-Box test
            GSPC_LB[j, k, l, m, i] <- GSPC_LBtest10$p.value
            
            GSPC_ARCH_LM[j, k, l, m, i] <- fDMA::archtest(GSPC_res, lag = 1)$p.value
            
            # Compute the Kolmogorov-Smirnov and Anderson-Darling statistics
            
            # Use the probability integral transform to transform to uniform random
            # variables
            
            # Use PIT to create uniform variables
            
            # For the 'norm' conditional distribution
            if (marginal_dists[i] == 'norm') {
              U <- pnorm(GSPC_res, mean=0, sd=1)
            }
            
            # For the 'snorm' conditional distribution
            if (marginal_dists[i] == 'snorm') {
              skew <- GSPC_model@fit$coef[['skew']]
              U <- fGarch::psnorm(GSPC_res, mean=0, sd=1, xi=skew)
            }
            
            # For the 'ged' conditional distribution
            if (marginal_dists[i] == 'ged') {
              shape <- GSPC_model@fit$coef[['shape']]
              U <- fGarch::pged(GSPC_res, mean=0, sd=1, nu=shape)
            }
            
            # For the 'sged' conditional distribution
            if (marginal_dists[i] == 'sged') {
              shape <- GSPC_model@fit$coef[['shape']]
              skew <- GSPC_model@fit$coef[['skew']]
              U <- fGarch::psged(GSPC_res, mean=0, sd=1, nu=shape, xi=skew)
            }
            
            # For the 'sstd' conditional distribution
            if (marginal_dists[i] == 'sstd') {
              shape <- GSPC_model@fit$coef[['shape']]
              skew <- GSPC_model@fit$coef[['skew']]
              U <- fGarch::psstd(GSPC_res, mean=0, sd=1, nu=shape, xi=skew)
            }
            
            # For the 'std' conditional distribution
            if (marginal_dists[i] == 'std') {
              shape <- GSPC_model@fit$coef[['shape']]
              U <- fGarch::pstd(GSPC_res, mean=0, sd=1, nu=shape)
            }
            
            # Conduct the KS test and record p-values
            GSPC_KS_AD[j, k, l, m, i, 1] <- KScorrect::LcKS(U, cdf='punif')$p.value
            
            # Conduct the AD test and record p-values
            GSPC_KS_AD[j, k, l, m, i, 2] <- ADGofTest::ad.test(U, null='punif')$p.value
          }
        }
      }
    }
  }
}

