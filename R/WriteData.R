library('quantmod')

##################################################
#
# Function to retrieve data and write it to the necessary files
#
##################################################

WriteLogRetPairs <- function(start, end, tickers, freq, file_path) {
  
  # Inputs
  #  - start : Start date for data retrieval in format %Y-%m-%d
  #  - end : End date for data retrieval in format %Y-%m-%d
  #  - tickers : Yahoo Finance tickers for which to retrieve data
  #  - freq : Frequency of the data, must be 'daily', 'weekly' or 'monthly'
  
  dataEnv <- new.env()
  
  start <- as.Date(start, format = '%Y-%m-%d')
  end <- as.Date(end, format = '%Y-%m-%d')
  
  getSymbols(
    Symbols = tickers,
    env = dataEnv,
    src = 'yahoo',
    from = start,
    to = end,
    periodicity = freq
  )
  
  write_tickers <- sort(ls(dataEnv), decreasing = F)
  n_tickers <- length(write_tickers)
  
  # Create data files following the format StartDate_EndDate_Index1_Index2_Freq_LogRets.csv
  for (i in 1:(n_tickers - 1)) {
    for (j in (i + 1):n_tickers) {
      file_name <- paste0(format.Date(start, format = '%Y%m%d'), '_',
                          format.Date(end, format = '%Y%m%d'), '_',
                          write_tickers[i], '_', write_tickers[j], '_', freq,
                          '_LogRets.csv')
      
      ticker1 <- get(write_tickers[i], envir = dataEnv)
      ticker2 <- get(write_tickers[j], envir = dataEnv)
      
      LogClose1 <- log(ticker1[ , paste0(write_tickers[i], '.Close')])
      LogClose2 <- log(ticker2[ , paste0(write_tickers[j], '.Close')])
      
      LogRets <- na.omit(diff(merge.xts(LogClose1, LogClose2, join = 'inner')))
      
      write.table(x = LogRets, file = paste0(file_path, file_name), row.names = F, col.names = F, sep = ',')
    }
  }
}

tickers <- c('^GSPC', '^FTSE', '^BVSP')
start <- '2016-01-01'
end <- '2022-05-31'
freq <- 'weekly'
file_path <- 'C://users/jonat/OneDrive/Documents/UCL MSc Statistics/STAT0034 MSc Project/Code/Data/'

WriteLogRetPairs(start, end, tickers, freq, file_path)

# Get dates for plotting

GetDates <- function(start, end, tickers, freq) {
  
  # Inputs
  #  - start : Start date for data retrieval in format %Y-%m-%d
  #  - end : End date for data retrieval in format %Y-%m-%d
  #  - tickers : Yahoo Finance tickers for which to retrieve data
  #  - freq : Frequency of the data, must be 'daily', 'weekly' or 'monthly'
  
  dataEnv <- new.env()
  
  start <- as.Date(start, format = '%Y-%m-%d')
  end <- as.Date(end, format = '%Y-%m-%d')
  
  quantmod::getSymbols(
    Symbols = tickers,
    env = dataEnv,
    src = 'yahoo',
    from = start,
    to = end,
    periodicity = freq
  )
  
  write_tickers <- sort(ls(dataEnv), decreasing = F)
  n_tickers <- length(write_tickers)
  
  # Find the common dates between the tickers after finding the log-returns
  Dates <- vector("list", n_tickers)
  for (i in 1:n_tickers) {
    ticker <- get(write_tickers[i], envir = dataEnv)
    LogClose <- log(ticker[ , paste0(write_tickers[i], '.Close')])
    Dates[[i]] <- zoo::index(na.omit(diff(LogClose)))
  }
  
  date_intersect <- Dates[[1]]
  for (i in 2:n_tickers) {
    date_intersect <- intersect(date_intersect, Dates[[i]])
  }
  
  return(as.Date(date_intersect, origin = "1970-01-01"))
}

