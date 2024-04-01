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

# Do the marginal modelling
n <- length(FTSE$FTSE.LogReturn)

GSPCLogRet <- GSPC$GSPC.LogReturn[2:n]
FTSELogRet <- FTSE$FTSE.LogReturn[2:n]

n <- length(FTSELogRet)

# Create plot of log-returns

library(RColorBrewer)

colours <- brewer.pal(4, 'Set1')

fig_dir <- 'C:/Users/jonat/OneDrive/Documents/UCL MSc Statistics/STAT0034 MSc Project/Report figures'

log_ret_fig_filename <- paste0(fig_dir, '/log_rets_fig.png')
png(log_ret_fig_filename, width = 1200, height = 1200, pointsize = 16)
par(mfrow = c(2, 1), cex = 1.15, cex.main = 1.5)
plot(as.numeric(GSPCLogRet) ~ zoo::index(GSPCLogRet), type = 'h', col = 'black',
     xlab = "Date", ylab = "Log-returns", bty = 'l',
     main = 'Daily log-returns for the S&P 500')
plot(as.numeric(FTSELogRet) ~ zoo::index(FTSELogRet), type = 'h', col = colours[2],
     xlab = 'Date', ylab = "Log-returns", bty = 'l',
     main = 'Daily log-returns for the FTSE 100')
dev.off()
  
  