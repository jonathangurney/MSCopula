# Set the working directory as desired
setwd('~/UCL MSc Statistics/STAT0034 MSc Project/Code')

# Source functions
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
source("./MSCopula/R/FitMSCopulaEM.R")

indices <- c('BVSP', 'FTSE', 'GSPC')
yahoo_indices <- c('^BVSP', '^FTSE', '^GSPC')
n_indices <- length(indices)

start_date <- as.Date("2000-01-01", format = "%Y-%m-%d")
end_date <- as.Date("2009-12-31", format = "%Y-%m-%d")
freq <- 'weekly'

# Set data directory
data_directory <- './Data'

# Set main figures directory, creating if necessary
figures_main_dir <- './Figures'
if (!(dir.exists(figures_main_dir))) {
  dir.create(figures_main_dir)
}

# Create directory for figures for this start_date, end_date combination
figures_dates_dir <- paste0(figures_main_dir, '/',
                                  format(start_date, format = "%Y%m%d"), '_', format(end_date, format = "%Y%m%d"))
if (!(dir.exists(figures_dates_dir))) {
  dir.create(figures_dates_dir)
}

# Set main logs directory, creating if necessary
log_dir <- './Logs'
if (!(dir.exists(log_dir))) {
  dir.create(log_dir)
}

optim_algo <- 'Nelder-Mead'
EM_max_iter <- 1e5 

for (i in 1:(n_indices - 1)) {
  for (j in (i + 1):n_indices) {
    # Set up directory for saving figures
    fig_dir <- paste0(figures_dates_dir, '/', indices[i], '_', indices[j], '_', freq)
    if (!(dir.exists(fig_dir))) {
      dir.create(fig_dir)
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
    
    par_init1 <- c(1, -1, 0, 1, 4, 5)
    theta0Gaussian <- c(0.5, 0.5)
    theta0Gumbel <- c(2, 2)
    theta0Clayton <- c(1, 1)
    theta0Frank <- c(1, 1)
    
    # Run Gaussian copula
    EMresultsGaussian <- FitMSCopulaEM(U, par_init1, theta0Gaussian, "gaussian", max_iter = EM_max_iter, algorithm = optim_algo)
    WriteFigs(EMresultsGaussian, fig_dir, indices[i], indices[j])
    WriteLogs(EMresultsGaussian, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)
    
    # Run Gumbel copula
    EMresultsGumbel <- FitMSCopulaEM(U, par_init1, theta0Gumbel, "gumbel", max_iter = EM_max_iter, algorithm = optim_algo)
    WriteFigs(EMresultsGumbel, fig_dir, indices[i], indices[j])
    WriteLogs(EMresultsGumbel, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)
    
    # Run Clayton copula
    EMresultsClayton <- FitMSCopulaEM(U, par_init1, theta0Clayton, "clayton", max_iter = EM_max_iter, algorithm = optim_algo)
    WriteFigs(EMresultsClayton, fig_dir, indices[i], indices[j])
    WriteLogs(EMresultsClayton, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)
    
    # Run Frank copula
    EMresultsFrank <- FitMSCopulaEM(U, par_init1, theta0Frank, "frank", max_iter = EM_max_iter, algorithm = optim_algo)
    WriteFigs(EMresultsFrank, fig_dir, indices[i], indices[j])
    WriteLogs(EMresultsFrank, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)
    
    par_init2 <- c(1, -1, 0, 1, 4, 5, 4, 4)
    theta0t <- c(0.5, 0.5)
    
    # Run t copula
    EMresultst <- FitMSCopulaEM(U, par_init2, theta0t, "t", max_iter = EM_max_iter, algorithm = optim_algo)
    WriteFigs(EMresultst, fig_dir, indices[i], indices[j])
    WriteLogs(EMresultst, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)
    
    par_init3 <- c(1, -1, 0, 1, 1, -1, 0, 1, 5, 4)
    theta0SJC <- c(0.5, 0.5)
    
    # Run SJC copula
    EMresultsSJC <- FitMSCopulaEM(U, par_init3, theta0SJC, "sjc", max_iter = EM_max_iter, algorithm = optim_algo)
    WriteFigs(EMresultsSJC, fig_dir, indices[i], indices[j])
    WriteLogs(EMresultsSJC, log_dir, start_date, end_date, indices = c(indices[i], indices[j]), freq = freq)
  }
}
