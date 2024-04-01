# Create simulated data and associated log files for the array jobs

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
source("./MSCopula/R/WriteJobLogs.R")

# Set the job id
job_id <- 1

# Get environment variables
env_vars <- Sys.getenv()
task_id <- env_vars[['SGE_TASK_ID']]

# Check if folder exists for this job, creating if necessary
job_dir <- paste0('./SimResults/', job_id)
if (!dir.exists(job_dir)) {
  dir.create(job_dir)
}

# Simulate a random number to ensure different seeds
start_time <- Sys.time()
rnorm_sim <- rnorm(1)

# Set some initial parameters
observations = 10
regimes = 2
family = "frank"

# Set list of parameters for simulating data
pars_list <- list(gaussian = c(2.5, 1, 0.3, 0.4),
                  gumbel = c(1.7, 0.7, 0.1, -1),
                  clayton = c(2, 1.4, 0.1, -1),
                  frank = c(3, 1.5, 0.3, -1.1))

# Set list of theta0 values
theta0_list <- list(gaussian = rep(0.5, regimes),
                    gumbel = rep(2, regimes),
                    clayton = rep(1, regimes),
                    frank = rep(1, regimes))

pars = pars_list[[family]]
pars_init = c(1, 2, 0, 1, 4, 5)

trans_mat = matrix(c(0.9, 0.1, 0.2, 0.8), ncol = 2)

theta0 = theta0_list[[family]]

# Create list of job settings
job_options <- list(job_id = job_id,
                    observations = observations,
                    regimes = regimes,
                    family = family,
                    pars = pars,
                    pars_init = pars_init,
                    trans_mat = trans_mat)

# Create log file for the job
WriteJobLogs(job_options, job_dir)

# Create seed file for this task
seed_dir <- paste0(job_dir, '/Seeds')
if (!dir.exists(seed_dir)) {
  dir.create(seed_dir)
}
seed_test <- matrix(c(format(start_time, "%Y-%m-%d %H:%M:%OS"), rnorm_sim), nrow = 1)
write.csv(seed_test, paste0(seed_dir, '/', job_id, '_', task_id, '_seed_test.csv'),
          row.names = F)

# Simulate the data
data <- SimData(n = observations, pars = pars,
                trans_mat = trans_mat, family = family)

# Create directory for storing the data
data_dir <- paste0(job_dir, '/Data')
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

# Write the data to csv files
write.table(data$U,
            paste0(data_dir, '/', job_id, '_', task_id, '_data.csv'),
            row.names = F, col.names = F, sep = ',')
write.table(data$states,
            paste0(data_dir, '/', job_id, '_', task_id, '_states.csv'),
            row.names = F, col.names = F, sep = ',')
write.table(data$theta,
            paste0(data_dir, '/', job_id, '_', task_id, '_theta.csv'),
            row.names = F, col.names = F, sep = ',')
write.table(data$theta_sim,
            paste0(data_dir, '/', job_id, '_', task_id, '_theta_sim.csv'),
            row.names = F, col.names = F, sep = ',')
