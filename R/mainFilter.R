# Set the working directory as desired
setwd('~/UCL MSc Statistics/STAT0034 MSc Project/Code')

# Source code
source("./MSCopula/R/gaussianCopulaDF.R")
source("./MSCopula/R/gumbelCopulaDF.R")
source("./MSCopula/R/claytonCopulaDF.R")
source("./MSCopula/R/frankCopulaDF.R")
source("./MSCopula/R/tCopulaDF.R")
source("./MSCopula/R/JCCopulaDF.R")
source("./MSCopula/R/SJCCopulaDF.R")
source("./MSCopula/R/pars2matrix.R")
source("./MSCopula/R/WriteData.R")
source("./MSCopula/R/WriteLogs.R")
source("./MSCopula/R/WriteFigs.R")
source("./MSCopula/R/MSCopulaLogLik.R")
source("./MSCopula/R/FitMSCopulaFilter.R")

# Set default optimization algorithm to be used
algorithm <- 'Nelder-Mead'

indices <- c('BVSP', 'FTSE', 'GSPC')
yahoo_indices <- c('^BVSP', '^FTSE', '^GSPC')
n_indices <- length(indices)

start_date <- as.Date("2000-01-01", format = "%Y-%m-%d")
end_date <- as.Date("2009-12-31", format = "%Y-%m-%d")
freq <- 'weekly'

# Set data directory
data_directory <- './Data'

# Set main figures directory, creating if necessary
figures_main_directory <- './Figures/Filter'
if (!(dir.exists(figures_main_directory))) {
  dir.create(figures_main_directory)
}

# Create directory for figures for this start_date, end_date combination
figures_dates_directory <- paste0(figures_main_directory, '/',
                                  format(start_date, format = "%Y%m%d"), '_', format(end_date, format = "%Y%m%d"), '_Figures/')
if (!(dir.exists(figures_dates_directory))) {
  dir.create(figures_dates_directory)
}

# Set main logs directory, creating if necessary
log_dir <- './Logs'
if (!(dir.exists(log_dir))) {
  dir.create(log_dir)
}

for (i in 1:(n_indices - 1)) {
  for (j in (i + 1):n_indices) {
    # Set up directory for saving figures
    fig_directory <- paste0(figures_dates_directory, indices[i], '_', indices[j], '_', freq, '/')
    if (!(dir.exists(fig_directory))) {
      dir.create(fig_directory)
    }
    
    data_filename <- paste0(format(start_date, format = "%Y%m%d"), '_', format(end_date, format = "%Y%m%d"), '_', indices[i], '_', indices[j], '_', freq, '_LogRets.csv')
    data <- read.csv(paste0(data_directory, '/', data_filename), header = FALSE)
    
    dates <- GetDates(start = format(start_date, format = '%Y-%m-%d'), format(end_date, format = '%Y-%m-%d'), tickers = yahoo_indices[c(i, j)], freq = freq)
    
    X <- as.matrix(data)
    empiricalCDF1 <- stats::ecdf(X[ , 1])
    empiricalCDF2 <- stats::ecdf(X[ , 2])
    
    U <- matrix(nrow = nrow(X), ncol = ncol(X))
    U[ , 1] <- 0.999 * empiricalCDF1(X[ , 1])
    U[ , 2] <- 0.999 * empiricalCDF2(X[ , 2])
    
    pars_init1 <- c(1, 2, 0, 1, 4, 5)
    theta0Gaussian <- c(0.5, 0.5)
    theta0Gumbel <- c(3, 3)
    theta0Clayton <- c(2, 2)
    theta0Frank <- c(1, 1)
    
    # Run Gaussian copula
    print(paste0('Fitting Gaussian copula for ', indices[i], ', ', indices[j]))
    GaussianResults <- FitMSCopulaFilter(U = U, pars_init = pars_init1, theta0 = theta0Gaussian, family = "gaussian", algorithm = algorithm)
    WriteFigs(GaussianResults, fig_directory, index1 = indices[i], index2 = indices[j])
    WriteLogs(GaussianResults, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)

    # Run Gumbel copula
    print(paste0('Fitting Gumbel copula for ', indices[i], ', ', indices[j]))
    GumbelResults <- FitMSCopulaFilter(U = U, pars_init = pars_init1, theta0 = theta0Gumbel, family = "gumbel", algorithm = algorithm)
    WriteFigs(GumbelResults, fig_directory, index1 = indices[i], index2 = indices[j])
    WriteLogs(GumbelResults, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)

    # Run Clayton copula
    print(paste0('Fitting Clayton copula for ', indices[i], ', ', indices[j]))
    ClaytonResults <- FitMSCopulaFilter(U = U, pars_init = pars_init1, theta0 = theta0Clayton, family = "clayton", algorithm = algorithm)
    WriteFigs(ClaytonResults, fig_directory, index1 = indices[i], index2 = indices[j])
    WriteLogs(ClaytonResults, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)
      
    # Run Frank copula
    print(paste0('Fitting Frank copula for ', indices[i], ', ', indices[j]))
    FrankResults <- FitMSCopulaFilter(U = U, pars_init = pars_init1, theta0 = theta0Frank, family = "frank", algorithm = algorithm)
    WriteFigs(FrankResults, fig_directory, index1 = indices[i], index2 = indices[j])
    WriteLogs(FrankResults, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)

    pars_init2 <- c(1, -1, 0, 1, 3, 3, 5, 4)
    theta0t <- c(0.5, 0.5)

    # Run t copula
    print(paste0('Fitting t copula for ', indices[i], ', ', indices[j]))
    tResults <- FitMSCopulaFilter(U = U, pars_init = pars_init2, theta0 = theta0t, family = "t", algorithm = algorithm)
    WriteFigs(tResults, fig_directory, index1 = indices[i], index2 = indices[j])
    WriteLogs(tResults, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)

    pars_init3 <- c(1, -1, 0, 1, 1, -1, 0, 1, 5, 4)
    theta0SJC <- c(0.5, 0.5)

    # Run SJC copula
    print(paste0('Fitting SJC copula for ', indices[i], ', ', indices[j]))
    SJCResults <- FitMSCopulaFilter(U = U, pars_init = pars_init3, theta0 = theta0SJC, family = "sjc", algorithm = algorithm)
    WriteFigs(SJCResults, fig_directory, index1 = indices[i], index2 = indices[j])
    WriteLogs(SJCResults, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)
  }
}
