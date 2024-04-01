setwd('~/Scratch/Code')
# setwd('~/UCL MSc Statistics/STAT0034 MSc Project/Code')

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
source("./MSCopula/R/FitMSCopulaEM.R")
source("./MSCopula/R/FitMSCopulaFilter.R")
source("./MSCopula/R/FitSimData.R")

# Set job id
job_id <- 4

print('Getting environment variables...')

# Get environment variables
env_vars <- Sys.getenv()
task_id <- env_vars[['SGE_TASK_ID']]

# Set job directory
job_dir <- paste0('./SimResults/', job_id)

print('Setting initial parameters...')

# Set some initial parameters
regimes = 2
family = "gumbel"

# Set list of theta0 values
theta0_list <- list(gaussian = rep(0.5, regimes),
                    gumbel = rep(2, regimes),
                    clayton = rep(1, regimes),
                    frank = rep(1, regimes))

theta0 = theta0_list[[family]]

pars_init = c(1, 2, 0, 1, 4, 5)

control <- list(HF_optim_algo = "Nelder-Mead",
                HF_NM_rel_tol = 1e-15,
                EM_optim_algo = "Nelder-Mead",
                EM_NM_rel_tol = 1e-9,
                EM_step_tol = 1e-3)

# Set data directory
data_dir <- paste0(job_dir, '/Data')

print('Reading data...')

# Read the data files
U <- as.matrix(read.csv(paste0(data_dir, '/', job_id, '_', task_id, '_data.csv'), header = F, sep = ','))
theta <- as.matrix(read.csv(paste0(data_dir, '/', job_id, '_', task_id, '_theta.csv'), header = F, sep = ','))
theta_sim <- as.matrix(read.csv(paste0(data_dir, '/', job_id, '_', task_id, '_theta_sim.csv'), header = F, sep = ','))
states <- as.matrix(read.csv(paste0(data_dir, '/', job_id, '_', task_id, '_states.csv'), header = F, sep = ','))

print("Fitting via Kim's filter...")

# Fit the MS Copula model using the HF algorithm
HFResult <- FitMSCopulaFilter(U = U,
                              pars_init = pars_init,
                              theta0 = theta0,
                              family = family,
                              regimes = regimes,
                              rel_tol = control$HF_NM_rel_tol)

# print("Fitting via EM...")
# 
# # Fit the MS copula model using the EM algorithm
# EMResult <- FitMSCopulaEM(U = U,
#                           pars_init = pars_init,
#                           theta0 = theta0,
#                           family = family,
#                           regimes = regimes,
#                           max_iter = 1e5,
#                           algorithm = control$EM_optim_algo,
#                           step_tol = control$EM_step_tol)

# Processing the data for output

n <- nrow(U)
d <- ncol(U)

print("Processing data...")

# Find the state with highest probability at each time period
HF_states <- apply(X = HFResult$eta, MARGIN = 1, FUN = function(x) {sample(which(x == max(x)), 1)})
# EM_states <- apply(X = EMResult$eta, MARGIN = 1, FUN = function(x) {sample(which(x == max(x)), 1)})

# Find the dependence parameter corresponding to the state with the highest
# probability at each time
HF_theta <- sapply(X = 1:n, FUN = function(x) {return(HFResult$theta[[1]][x, HF_states[x]])})
# EM_theta <- sapply(X = 1:n, FUN = function(x) {return(EMResult$theta[[1]][x, EM_states[x]])})

# Find the mean error, RMSE, MAE
HF_ME <- mean(HF_theta - theta_sim)
# EM_ME <- mean(EM_theta - theta_sim)
HF_RMSE <- sqrt(mean((HF_theta - theta_sim)^2))
# EM_RMSE <- sqrt(mean((EM_theta - theta_sim)^2))
HF_MAE <- mean(abs(HF_theta - theta_sim))
# EM_MAE <- mean(abs(EM_theta - theta_sim))

# Find the log-likelihood, AIC, and BIC for the fitted models
HF_loglik <- HFResult$loglik
# EM_loglik <- EMResult$loglik
HF_AIC <- HFResult$aic
# EM_AIC <- EMResult$aic
HF_BIC <- HFResult$bic
# EM_BIC <- EMResult$bic

# Create csv filenames
HF_errors_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, '_HF_errors_.csv')
# EM_errors_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, '_EM_errors_.csv')
HF_pars_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, '_HF_parameters_.csv')
# EM_pars_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, '_EM_parameters_.csv')
HF_GOF_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, '_HF_GOF_.csv')
# EM_GOF_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, '_EM_GOF_.csv')
runtimes_filename <- paste0(data_dir, '/', job_id, '_', task_id, '_', family, '_runtimes_.csv')

filenames <- c(HF_errors_filename, HF_pars_filename, HF_GOF_filename, runtimes_filename)

print('Creating dataframes...')

HF_errors <- matrix(c(HF_ME, HF_RMSE, HF_MAE), nrow = 1)
# EM_errors <- matrix(c(EM_ME, EM_RMSE, EM_MAE), nrow = 1)
HF_pars <- matrix(HFResult$pars, nrow = 1)
# EM_pars <- matrix(EMResult$pars, nrow = 1)
HF_GOF <- matrix(c(HF_loglik, HF_AIC, HF_BIC), nrow = 1)
# EM_GOF <- matrix(c(EM_loglik, EM_AIC, EM_BIC), nrow = 1)
runtimes <- matrix(c(HFResult$runtime), nrow = 1)

# Write output dataframes to files

print('Writing output...')

write.table(HF_errors, file = HF_errors_filename, row.names = F, col.names = F, sep = ',')
# write.table(EM_errors, file = EM_errors_filename, row.names = F, col.names = F, sep = ',')
write.table(HF_pars, file = HF_pars_filename, row.names = F, col.names = F, sep = ',')
# write.table(EM_pars, file = EM_pars_filename, row.names = F, col.names = F, sep = ',')
write.table(HF_GOF, file = HF_GOF_filename, row.names = F, col.names = F, sep = ',')
# write.table(EM_GOF, file = EM_GOF_filename, row.names = F, col.names = F, sep = ',')
write.table(runtimes, file = runtimes_filename, row.names = F, col.names = F, sep = ',')

