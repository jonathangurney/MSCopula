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
  common_dates <- as.Date(intersect(zoo::index(FTSE), zoo::index(GSPC)), origin = "1970-01-01")
  GSPC <- GSPC[common_dates, ]
  FTSE <- FTSE[common_dates, ]
}

# Create vector of log close prices
GSPC$GSPC.LogClose <- log(GSPC$GSPC.Close)
FTSE$FTSE.LogClose <- log(FTSE$FTSE.Close)

# Create vector of log returns
GSPC$GSPC.LogReturn <- c(diff(GSPC$GSPC.LogClose))
FTSE$FTSE.LogReturn <- c(diff(FTSE$FTSE.LogClose))

n <- length(FTSE$FTSE.LogReturn)

GSPCLogRet <- GSPC$GSPC.LogReturn[2:n]
FTSELogRet <- FTSE$FTSE.LogReturn[2:n]

# Plot the acf and pacf for the returns and the squared returns
max_lag <- 25

fig_dir <- paste0('C:/Users/jonat/OneDrive/Documents/UCL MSc Statistics/STAT0034 MSc Project/Report figures')

FTSE_ret_acf_filename <- paste0(fig_dir, '/FTSE_rets_sq_rets_acf.png')
GSPC_ret_acf_filename <- paste0(fig_dir, '/GSPC_rets_sq_rets_acf.png')
FTSE_sq_ret_acf_filename <- paste0(fig_dir, '/FTSE_sq_rets_acf.png')
GSPC_sq_ret_acf_filename <- paste0(fig_dir, '/GSPC_sq_rets_acf.png')

png(filename = FTSE_ret_acf_filename, width = 1200, height = 800, pointsize = 16)
par(mfrow = c(2, 2), cex = 1.15)
acf(FTSELogRet, col = 'green', lwd = 2, lag.max = max_lag,
    main = 'ACF for log-returns of FTSE 100 index')
pacf(FTSELogRet, col = 'red', lwd = 2, lag.max = max_lag,
     main = 'PACF for log-returns of FTSE 100 index')
acf(FTSELogRet^2, col = 'green', lwd = 2, lag.max = max_lag,
    main = 'ACF for squared log-returns of FTSE 100 index')
pacf(FTSELogRet^2, col = 'red', lwd = 2, lag.max = max_lag,
     main = 'PACF for squared log-returns of FTSE 100 index')
dev.off()

png(filename = GSPC_ret_acf_filename, width = 1200, height = 800, pointsize = 16)
par(mfrow = c(2, 2), cex = 1.15)
acf(GSPCLogRet, col = 'green', lwd = 2, lag.max = max_lag,
    main = 'ACF for log-returns of S&P 500 index')
pacf(GSPCLogRet, col = 'red', lwd = 2, lag.max = max_lag,
     main = 'PACF for log-returns of S&P 500 index')
acf(GSPCLogRet^2, col = 'green', lwd = 2, lag.max = max_lag,
    main = 'ACF for squared log-returns of S&P 500 index')
pacf(GSPCLogRet^2, col = 'red', lwd = 2, lag.max = max_lag,
     main = 'PACF for squared log-returns of S&P 500 index')
dev.off()

# Create tables showing parameter estimates



