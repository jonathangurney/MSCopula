# setwd('~/Scratch/Code')
setwd('~/UCL MSc Statistics/STAT0034 MSc Project/Code')

source("./MSCopula/R/gaussianCopulaDF.R")
source("./MSCopula/R/gumbelCopulaDF.R")
source("./MSCopula/R/claytonCopulaDF.R")
source("./MSCopula/R/frankCopulaDF.R")
source("./MSCopula/R/tCopulaDF.R")
source("./MSCopula/R/JCCopulaDF.R")
source("./MSCopula/R/SJCCopulaDF.R")
source("./MSCopula/R/pars2matrix.R")
source("./MSCopula/R/WriteLogs.R")
source("./MSCopula/R/WriteFigs.R")
source("./MSCopula/R/SimData.R")
source("./MSCopula/R/MSCopulaLogLik.R")
# source("./MSCopula/R/FitMSCopulaEM.R")
source("./MSCopula/R/FitMSCopulaFilter.R")
source("./MSCopula/R/FitSimData.R")

# Set job id
job_id <- 10

print('Getting environment variables...')

# Get environment variables
env_vars <- Sys.getenv()
task_id <- 1

# Set job directory
job_dir <- paste0('./SimResults/', job_id)
if (!dir.exists(job_dir)) {
  dir.create(job_dir)
}

print('Setting initial parameters...')

wallclock_time <- 1.1

# Do we get initial parameters from a checkpoint file?
checkpoint <- T
checkpoint_dir <- paste0(job_dir, '/Checkpoint')
if (!dir.exists(checkpoint_dir)) {
  stop('Checkpoint is TRUE but directory does not exist.')
}

checkpoint_filename <- paste0(checkpoint_dir, '/', job_id, '_', task_id, '_', family, '_checkpoint_data.csv')
checkpoint_data <- as.matrix(read.csv(checkpoint_filename, header = F))

# Set some initial parameters
regimes = 2
family = "gumbel"

# Set list of theta0 values
theta0_list <- list(gaussian = rep(0.5, regimes),
                    gumbel = rep(2, regimes),
                    clayton = rep(1, regimes),
                    frank = rep(1, regimes))

theta0 = theta0_list[[family]]

if (checkpoint) {
  pars_init <- as.numeric(checkpoint_data[1:6])
  init_runtime <- checkpoint_data[7]
} else {
  pars_init = c(1, 2, 0, 1, 4, 5)
  init_runtime <- 0
}

control <- list(EM_optim_algo = "Nelder-Mead",
                EM_NM_rel_tol = 1e-9,
                EM_step_tol = 1e-1)

# Set data directory
data_dir <- paste0(job_dir, '/Data')
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

print('Reading data...')

set.seed(1)

# Read the data files
U <- matrix(runif(20), ncol = 2)

# theta <- as.matrix(read.csv(paste0(data_dir, '/', job_id, '_', task_id, '_theta.csv'), header = F, sep = ','))
# theta_sim <- as.matrix(read.csv(paste0(data_dir, '/', job_id, '_', task_id, '_theta_sim.csv'), header = F, sep = ','))
# states <- as.matrix(read.csv(paste0(data_dir, '/', job_id, '_', task_id, '_states.csv'), header = F, sep = ','))


print("Fitting via EM...")

# Fit the MS copula model using the EM algorithm
EMResult <- FitMSCopulaEM(U = U,
                          pars_init = pars_init,
                          theta0 = theta0,
                          family = family,
                          regimes = regimes,
                          max_iter = 1e5,
                          algorithm = control$EM_optim_algo,
                          step_tol = control$EM_step_tol,
                          wall_time = wallclock_time)

if (!EMResult$convergence) {
  # Write parameter values to file for checkpointing
  EM_checkpoint <- matrix(c(EMResult$pars, (init_runtime + EMResult$runtime)), nrow = 1)

  if (!dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir)
  }
  
  write.table(EM_checkpoint, checkpoint_filename, row.names = F, col.names = F, sep = ',')
  
} else {
  # Processing the data for output
  
  n <- nrow(U)
  d <- ncol(U)
  
  print("Processing data...")
  
  # Find the state with highest probability at each time period
  EM_states <- apply(X = EMResult$eta, MARGIN = 1, FUN = function(x) {sample(which(x == max(x)), 1)})
  
  # Find the dependence parameter corresponding to the state with the highest
  # probability at each time
  EM_theta <- sapply(X = 1:n, FUN = function(x) {return(EMResult$theta[[1]][x, EM_states[x]])})
  
  # Find the mean error, RMSE, MAE
  EM_ME <- mean(EM_theta - theta_sim)
  EM_RMSE <- sqrt(mean((EM_theta - theta_sim)^2))
  EM_MAE <- mean(abs(EM_theta - theta_sim))
  
  # Find the log-likelihood, AIC, and BIC for the fitted models
  EM_loglik <- EMResult$loglik
  EM_AIC <- EMResult$aic
  EM_BIC <- EMResult$bic
  
  # Create csv filenames
  EM_errors_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, '_EM_errors_.csv')
  EM_pars_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, '_EM_parameters_.csv')
  EM_GOF_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, '_EM_GOF_.csv')
  EM_runtimes_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, 'EM_runtimes_.csv')
  
  filenames <- c(EM_errors_filename, EM_pars_filename, EM_GOF_filename, EM_runtimes_filename)
  
  print('Creating dataframes...')
  
  EM_errors <- matrix(c(EM_ME, EM_RMSE, EM_MAE), nrow = 1)
  EM_pars <- matrix(EMResult$pars, nrow = 1)
  EM_GOF <- matrix(c(EM_loglik, EM_AIC, EM_BIC), nrow = 1)
  EM_runtimes <- matrix(c(init_runtime + EMResult$runtime), nrow = 1)
  
  # Write output dataframes to files
  
  print('Writing output...')
  
  write.table(EM_errors, file = EM_errors_filename, row.names = F, col.names = F, sep = ',')
  write.table(EM_pars, file = EM_pars_filename, row.names = F, col.names = F, sep = ',')
  write.table(EM_GOF, file = EM_GOF_filename, row.names = F, col.names = F, sep = ',')
  write.table(EM_runtimes, file = EM_runtimes_filename, row.names = F, col.names = F, sep = ',')
}
