setwd('~/Scratch/Code')

library(snow)

# Create reference to cluster
cl <- snow::getMPIcluster()

source("./MSCopula/R/FitSimData.R")

# Source code on each machine
snow::clusterEvalQ(cl, {
  library(parallel)
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
  source("./MSCopula/R/SimData.R")
  source("./MSCopula/R/MSCopulaLogLik.R")
  source("./MSCopula/R/FitMSCopulaEM.R")
  source("./MSCopula/R/FitMSCopulaFilter.R")
  source("./MSCopula/R/FitSimData.R")
})

# Set up RNG stream for each node
snow::clusterSetupRNG(cl)

# Initialize parameters
pars_list <- list(gaussian = c(2.5, 1, 0.3, 0.4),
                  gumbel = c(1.7, 0.7, 0.1, -1),
                  clayton = c(2, 1.4, 0.1, -1),
                  frank = c(3, 1.5, 0.3, -1.1))

replications <- 10
observations <- 10

pars_init <- c(1, 2, 0, 1, 4, 5)
P <- matrix(c(0.99, 0.01, 0.025, 0.975), ncol = 2)
theta0 <- c(0.5, 0.5)
family <- "gaussian"

output <- clusterApplyLB(cl = cl, rep(observations, replications), fun = FitSimData, pars = pars_list[[family]],
               pars_init = pars_init, trans_mat = P, theta0 = theta0, family = family)

# Create dataframes for writing to files
HF_errors <- matrix(ncol = 3, nrow = replications)
EM_errors <- matrix(ncol = 3, nrow = replications)
HF_pars <- matrix(ncol = length(pars_init), nrow = replications)
EM_pars <- matrix(ncol = length(pars_init), nrow = replications)
runtimes <- nmatrix(nrow = replications, ncol = 2)

for (i in 1:replications) {
  errors[i, ] <- cbind(output[[i]][['HF_errors']], output[[i]][['EM_errors']])
  pars[i, ] <- cbind(output[[i]][['HF_pars']], output[[i]][['EM_pars']])
  runtimes[i, ] <- output[[i]][['runtimes']]
}

# Write output dataframes to files
errors_filename <- paste0('./Logs/', family, '_errors_', observations, '.csv')
pars_filename <- paste0('./Logs/', family, '_parameters_', observations, '.csv')
runtimes_filename <- paste0('./Logs/', family, '_runtimes_', observations, '.csv')
write.csv(errors, file = errors_filename, row.names = F)
write.csv(pars, file = pars_filename, row.names = F)
write.csv(runtimes, file = runtimes_filename, row.names = F)

snow::stopImplicitCluster()



