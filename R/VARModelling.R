library(quantmod)
library(stats)

start_date <- as.Date('01-01-1990', format = "%d-%m-%Y")
end_date <- as.Date('30-06-2022', format = "%d-%m-%Y")

getSymbols(Symbols = c('^FTSE', '^BVSP', '^GSPC'),
           src = 'yahoo',
           from = start_date,
           to = end_date,
           periodicity = 'weekly')

# Check if data contains any NA values

if (any(is.na(GSPC))) {
  stop('GSPC contains NA values')
}

if (any(is.na(FTSE))) {
  stop('FTSE contains NA values')
}

if (length(setdiff(zoo::index(FTSE), zoo::index(GSPC))) != 0) {
  common_dates <- intersect(zoo:index(FTSE), zoo::index(GSPC))
  GSPC <- GSPC[common_dates]
  FTSE <- FTSE[common_dates]
}

# Create vector of log close prices
GSPC$GSPC.LogClose <- log(GSPC$GSPC.Close)
FTSE$FTSE.LogClose <- log(FTSE$FTSE.Close)

# Create vector of log returns
GSPC$GSPC.LogReturn <- c(diff(GSPC$GSPC.LogClose))
FTSE$FTSE.LogReturn <- c(diff(FTSE$FTSE.LogClose))

# Do the marginal modelling

GSPCLogRet <- GSPC$GSPC.LogReturn[2:n]
FTSELogRet <- FTSE$FTSE.LogReturn[2:n]

n <- nrow(FTSELogRet)

# Plot the acf and pacf for the returns and the squared returns
max_lag <- 30
par(mfrow = c(2, 1))
acf(GSPCLogRet, col = 'green', lwd = 2, lag.max = max_lag)
pacf(GSPCLogRet, col = 'red', lwd = 2, lag.max = max_lag)
acf(FTSELogRet, col = 'green', lwd = 2, lag.max = max_lag)
pacf(FTSELogRet, col = 'red', lwd = 2, lag.max = max_lag)
acf(GSPCLogRet^2, col = 'green', lwd = 2, lag.max = max_lag)
pacf(GSPCLogRet^2, col = 'red', lwd = 2, lag.max = max_lag)
acf(FTSELogRet^2, col = 'green', lwd = 2, lag.max = max_lag)
pacf(FTSELogRet^2, col = 'red', lwd = 2, lag.max = max_lag)

# Models to consider for the FTSE
FTSE_ARMA_orders <- list(ar = c(1, 3, 7), ma = c(1, 3, 7))
FTSE_GARCH_orders <- list(ar = c(1, 2, 3), ma = c(1, 2, 3))

GSPC_ARMA_orders <- list(ar = c(1), ma = c(1))
GSPC_GARCH_orders <- list(ar = c(1, 2, 3), ma = c(1, 2, 3))

# Get indices
n_FTSE_p <- length(FTSE_ARMA_orders$ar)
n_FTSE_q <- length(FTSE_ARMA_orders$ma)
n_FTSE_r <- length(FTSE_GARCH_orders$ar)
n_FTSE_s <- length(FTSE_GARCH_orders$ma)

# Initialize arrays for storing GOF test results

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
FTSE_LB <- array(dim = c(FTSE_array_dims, 4))

# Initialize array for storing Kolmogorov-Smirnov and Anderson-Darling p-values
FTSE_KS_AD <- array(dim = c(FTSE_array_dims, 2))

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
            FTSE_IC[j, k, l, m, i, ] <- c(FTSE_model@fit$ics[['AIC']], FTSE_model@fit$ics[['BIC']])

            # Compute Ljung-Box p-value for residuals and squared residuals
            FTSE_LBtest10 <- stats::Box.test(FTSE_res, lag = 10, type = 'Ljung-Box', fitdf = (p + q + r + s))
            FTSE_LBtest10_2 <- stats::Box.test(FTSE_res^2, lag = 10, type = 'Ljung-Box', fitdf = (p + q + r + s))
            FTSE_LBtest20 <- stats::Box.test(FTSE_res, lag = 20, type = 'Ljung-Box', fitdf = (p + q + r + s))
            FTSE_LBtest20_2 <- stats::Box.test(FTSE_res^2, lag = 20, type = 'Ljung-Box', fitdf = (p + q + r + s))

            #Record p-values fromt the Ljung-Box test
            FTSE_LB[j, k, l, m, i, ] <- c(FTSE_LBtest10$p.value,
                                          FTSE_LBtest10_2$p.value,
                                          FTSE_LBtest20$p.value,
                                          FTSE_LBtest20_2$p.value)

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

FTSE_LB1 <- FTSE_LB[ , , , , , 1]
FTSE_LB2 <- FTSE_LB[ , , , , , 2]
FTSE_LB3 <- FTSE_LB[ , , , , , 3]
FTSE_LB4 <- FTSE_LB[ , , , , , 4]

FTSE_KS <- FTSE_KS_AD[ , , , , , 1]
FTSE_AD <- FTSE_KS_AD[ , , , , , 2]

plot(FTSE_KS)
plot(FTSE_AD)
colour <- c("red", "green")[(FTSE_LB1[ , , , , 4:6] > 0.05) + 1]
plot(FTSE_AD[ , , , , 4:6], col = colour)

#############################################################################
#
# Same but now it's GSPC
#
#############################################################################

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
GSPC_LB <- array(dim = c(GSPC_array_dims, 4))

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
            GSPC_IC[j, k, l, m, i, ] <- c(GSPC_model@fit$ics[['AIC']], GSPC_model@fit$ics[['BIC']])
            
            # Compute Ljung-Box p-value for residuals and squared residuals
            GSPC_LBtest10 <- stats::Box.test(GSPC_res, lag = 10, type = 'Ljung-Box', fitdf = (p + q + r + s))
            GSPC_LBtest10_2 <- stats::Box.test(GSPC_res^2, lag = 10, type = 'Ljung-Box', fitdf = (p + q + r + s))
            GSPC_LBtest20 <- stats::Box.test(GSPC_res, lag = 20, type = 'Ljung-Box', fitdf = (p + q + r + s))
            GSPC_LBtest20_2 <- stats::Box.test(GSPC_res^2, lag = 20, type = 'Ljung-Box', fitdf = (p + q + r + s))
            
            #Record p-values fromt the Ljung-Box test
            GSPC_LB[j, k, l, m, i, ] <- c(GSPC_LBtest10$p.value,
                                          GSPC_LBtest10_2$p.value,
                                          GSPC_LBtest20$p.value,
                                          GSPC_LBtest20_2$p.value)
            
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

GSPC_LB1 <- GSPC_LB[ , , , , , 1]
GSPC_LB2 <- GSPC_LB[ , , , , , 2]
GSPC_LB3 <- GSPC_LB[ , , , , , 3]
GSPC_LB4 <- GSPC_LB[ , , , , , 4]

GSPC_KS <- GSPC_KS_AD[ , , , , , 1]
GSPC_AD <- GSPC_KS_AD[ , , , , , 2]

plot(GSPC_KS)
plot(GSPC_AD)
colour <- c("red", "green")[(GSPC_LB1[ , , 1:6] > 0.05) + 1]
plot(GSPC_AD[ , , 1:6], col = colour)

