setwd('~/Scratch/Code')

library(snow)

print("Getting cluster...")

# Create reference to cluster
cl <- snow::getMPIcluster()

print("Got cluster")
print("Sourcing modules on main")

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

print("Sourced modules on main")
print("Sourcing modules on cluster")

# Set wd of the cluster nodes
snow::clusterEvalQ(cl, setwd('~/Scratch/Code'))

# Source code on each machine
snow::clusterEvalQ(cl, source("./MSCopula/R/gaussianCopulaDF.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/gumbelCopulaDF.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/claytonCopulaDF.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/frankCopulaDF.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/tCopulaDF.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/JCCopulaDF.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/SJCCopulaDF.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/pars2matrix.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/WriteLogs.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/WriteFigs.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/SimData.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/MSCopulaLogLik.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/FitMSCopulaEM.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/FitMSCopulaFilter.R"))
snow::clusterEvalQ(cl, source("./MSCopula/R/FitSimData.R"))

print("Sourced modules on cluster")
print("Setting up RNG...")

# Set up RNG stream for each node
snow::clusterSetupRNG(cl)

print("Set up RNG")
print("Initializing parameters")

# Initialize parameters
pars_list <- list(gaussian = c(2.5, 1, 0.3, 0.4),
                  gumbel = c(1.7, 0.7, 0.1, -1),
                  clayton = c(2, 1.4, 0.1, -1),
                  frank = c(3, 1.5, 0.3, -1.1))

replications <- 100
observations <- 100

pars_init <- c(1, 2, 0, 1, 4, 5)
P <- matrix(c(0.99, 0.01, 0.025, 0.975), ncol = 2)
theta0 <- c(0.5, 0.5)
family <- "gaussian"
regimes <- 2

print("Running clusterApplyLB")

output <- clusterApplyLB(cl = cl, rep(observations, replications), fun = FitSimData, pars = pars_list[[family]], regimes = regimes,
                         pars_init = pars_init, trans_mat = P, theta0 = theta0, family = family)

print("Processing data")

# Create dataframes for writing to files
errors <- matrix(ncol = 6, nrow = replications)
pars <- matrix(ncol = (2 * length(pars_init)), nrow = replications)
runtimes <- matrix(nrow = replications, ncol = 2)

for (i in 1:replications) {
  errors[i, ] <- cbind(output[[i]][['HF_errors']], output[[i]][['EM_errors']])
  pars[i, ] <- cbind(output[[i]][['HF_pars']], output[[i]][['EM_pars']])
  runtimes[i, ] <- output[[i]][['runtimes']]
}

colnames(errors) <- c('HF_ME', 'HF_RMSE', 'HF_MAE', 'EM_ME', 'EM_RMSE', 'EM_MAE')
colnames(pars) <- c(paste0('HF_omega', 1:regimes), 'HF_alpha', 'HF_beta', 'HF_p11', 'HF_p22',
                    paste0('EM_omega', 1:regimes), 'EM_alpha', 'EM_beta', 'EM_p11', 'EM_p22')
colnames(runtimes) <- c('HF.runtime', 'EM.runtime')

print("Writing data...")

# Write output dataframes to files
errors_filename <- paste0('./SimResults/', family, '_errors_', observations, '.csv')
pars_filename <- paste0('./SimResults/', family, '_parameters_', observations, '.csv')
runtimes_filename <- paste0('./SimResults/', family, '_runtimes_', observations, '.csv')
write.csv(errors, file = errors_filename, row.names = F)
write.csv(pars, file = pars_filename, row.names = F)
write.csv(runtimes, file = runtimes_filename, row.names = F)

print("Data writing complete")

snow::stopCluster(cl)

