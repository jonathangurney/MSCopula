# setwd('~/Scratch/Code')
setwd('~/UCL MSc Statistics/STAT0034 MSc Project/Code')

# Load the necessary packages
library('fGarch')
library('ADGofTest')
library('KScorrect')

# Source the files to be able to fit the copula models to the data
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

# Set some initial parameters
regimes = 2
family = "t"

# Set the length of the rolling window
window_length <- 1000
observations <- window_length

# Set list of theta0 values
theta0_list <- list(gaussian = rep(0.5, regimes),
                    gumbel = rep(2, regimes),
                    clayton = rep(1, regimes),
                    frank = rep(1, regimes),
                    t = rep(0.5, regimes))

theta0 = theta0_list[[family]]

pars_init = c(1, 2, 0, 1, 4, 5, 4, 4)

control <- list(HF_optim_algo = "Nelder-Mead",
                HF_NM_rel_tol = 1e-12,
                EM_optim_algo = "Nelder-Mead",
                EM_NM_rel_tol = 1e-12,
                EM_step_tol = 2e-3)

port_weights <- c(0.5, 0.5)

n_MC <- 10000

# Get environment variables
env_vars <- Sys.getenv()
task_id <- as.numeric(1)

# Import the data from the VaRData directory
var_dir <- paste0('./VaRBacktest/Data')
var_data_filename <- paste0(var_dir, '/data.csv')
data <- read.csv(var_data_filename, header = F)
data <- data[task_id:(task_id + window_length - 1), ]

# Fit the AR-GARCH models using fGarch
FTSE_garch <- fGarch::garchFit(formula = ~ garch(1, 1),
                               data = data[ , 1],
                               cond.dist = 'sged',
                               trace = F)

GSPC_garch <- fGarch::garchFit(formula = ~ arma(1, 0) + garch(1, 1),
                               data = data[ , 2],
                               cond.dist = 'sged',
                               trace = F)

# Create vectors of residuals
FTSE_res <- fGarch::residuals(FTSE_garch, standardize = T)
GSPC_res <- fGarch::residuals(GSPC_garch, standardize = T)

# Use the probability integral transform to create uniform random variables
u1 <- fGarch::psged(FTSE_res,
                    xi = FTSE_garch@fit$coef['skew'],
                    nu = FTSE_garch@fit$coef['shape'])

u2 <- fGarch::psged(GSPC_res,
                    xi = GSPC_garch@fit$coef['skew'],
                    nu = GSPC_garch@fit$coef['shape'])

U <- cbind(u1, u2)

print(paste0('u1 KS test p-value = ', KScorrect::LcKS(u1, cdf='punif')$p.value))
print(paste0('u1 AD test p-value = ', ADGofTest::ad.test(u1, null = 'punif')$p.value))
print(paste0('u2 KS test p-value = ', KScorrect::LcKS(u2, cdf='punif')$p.value))
print(paste0('u2 AD test p-value = ', ADGofTest::ad.test(u2, null = 'punif')$p.value))

# Now fit the regime switching copula model using both IFM and EM

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

# Simplify notation with the parameters
HF_omega <- HFResult$pars[1:2]
HF_alpha <- HFResult$pars[3]
HF_beta <- HFResult$pars[4]
HF_trans_mat <- pars2matrix(HFResult$pars[5:6])

# EM_omega <- EMResult$pars[1:2]
# EM_alpha <- EMResult$pars[3]
# EM_beta <- EMResult$pars[4]
# EM_trans_mat <- pars2matrix(EMResult$pars[5:6])

# Compute the predictive probability
HF_pred_prob <- HFResult$eta[observations, ] %*% HF_trans_mat
# EM_pred_prob <- EMResult$eta[observations, ] %*% EM_trans_mat

# Compute the next value of the dependence parameter
switch(family, 
       "clayton" = {
         psi <- mean(abs(U[(observations - 9):observations, 1] - U[(observations - 9):observations, 2]))
         
         transform_func <- function(x) {
           return(x^2)
         }
       },
       
       "gumbel" = {
         psi <- mean(abs(U[(observations - 9):observations, 1] - U[(observations - 9):observations, 2]))
         
         transform_func <- function(x) {
           return(1 + x^2)
         }
       },
       
       "gaussian" = {
         psi <- mean(qnorm(U[(observations - 9):observations, 1]) * qnorm(U[(observations - 9):observations, 2]))
         
         transform_func <- function(x) {
           return(1.9998/(exp(-x) + 1) - 0.9999)
         }
       })

# Compute the parameter value at the next stage for both regimes
HF_lp <- HF_omega + HF_alpha * HFResult$theta[[1]][observations, ] + HF_beta * psi
# EM_lp <- EM_omega + EM_alpha * EMResult$theta[[1]][observations, ] + EM_beta * psi

HF_theta_new <- transform_func(HF_lp)
# EM_theta_new <- transform_func(EM_lp)

# Choose the predicted regimes and select the appropriate value of theta for
# each estimation method
HF_state <- sample(which(HF_pred_prob >= 0.5), size = 1)
# EM_state <- sample(which(EM_pred_prob >= 0.5), size = 1)

HF_theta_new <- HF_theta_new[HF_state]
# EM_theta_new <- EM_theta_new[EM_state]

# Simulate the equally weighted portfolio returns and compute 95% and 99% VaR
HF_uniform_sims <- matrix(nrow = n_MC, ncol = 2)
# EM_uniform_sims <- matrix(nrow = n_MC, ncol = 2)

# Create copulas with the relevant parameters for simulation
switch(family, 
       "clayton" = {
         HF_cop <- copula::claytonCopula(param = HF_theta_new)
         # EM_cop <- copula::claytonCopula(param = EM_theta_new)
       },
       
       "gumbel" = {
         HF_cop <- copula::gumbelCopula(param = HF_theta_new)
         # EM_cop <- copula::gumbelCopula(param = EM_theta_new)
       },
       
       "gaussian" = {
         HF_cop <- copula::normalCopula(param = HF_theta_new)
         # EM_cop <- copula::normalCopula(param = EM_theta_new)
       })

# Simulate values from the copulas
HF_uniform_sims <- copula::rCopula(n = n_MC, copula = HF_cop)
# EM_uniform_sims <- copula::rCopula(n = n_MC, copula = EM_cop)

# Use inverse probability integral transform to generate innovations
HF_FTSE_innov_sims <- fGarch::qsged(HF_uniform_sims[ , 1],
                                  nu = FTSE_garch@fit$coef['shape'],
                                  xi = FTSE_garch@fit$coef['skew'])

HF_GSPC_innov_sims <- fGarch::qsged(HF_uniform_sims[ , 2],
                                  nu = GSPC_garch@fit$coef['shape'],
                                  xi = GSPC_garch@fit$coef['skew'])

# EM_FTSE_innov_sims <- fGarch::qsged(EM_uniform_sims[ , 1],
#                                   nu = FTSE_garch@fit$coef['shape'],
#                                   xi = FTSE_garch@fit$coef['skew'])
# 
# EM_GSPC_innov_sims <- fGarch::qsged(EM_uniform_sims[ , 2],
#                                   nu = GSPC_garch@fit$coef['shape'],
#                                   xi = GSPC_garch@fit$coef['skew'])

# Compute the conditional mean and variance for the prediction time period
FTSE_mu <- FTSE_garch@fit$coef['mu']
FTSE_sigma <- sqrt(FTSE_garch@fit$coef['omega'] + FTSE_garch@fit$coef['alpha1'] * FTSE_garch@residuals[observations]^2 +
                 FTSE_garch@fit$coef['beta1'] * FTSE_garch@h.t[observations])

GSPC_mu <- GSPC_garch@fit$coef['mu'] + GSPC_garch@fit$coef['ar1'] * GSPC_garch@fitted[observations]
GSPC_sigma <- sqrt(GSPC_garch@fit$coef['omega'] + GSPC_garch@fit$coef['alpha1'] * GSPC_garch@residuals[observations]^2 +
                        GSPC_garch@fit$coef['beta1'] * GSPC_garch@h.t[observations])

# Compute the simulated marginal returns
HF_FTSE_ret_sims <- FTSE_mu + FTSE_sigma * HF_FTSE_innov_sims
HF_GSPC_ret_sims <- GSPC_mu + GSPC_sigma * HF_GSPC_innov_sims

# EM_FTSE_ret_sims <- FTSE_mu + FTSE_sigma * EM_FTSE_innov_sims
# EM_GSPC_ret_sims <- GSPC_mu + GSPC_sigma * EM_GSPC_innov_sims

# Compute simulated portfolio returns and 95% and 99% portfolio VaR
HF_port_ret_sims <- (exp(cbind(HF_FTSE_ret_sims, HF_GSPC_ret_sims)) %*% port_weights) - 1
# EM_port_ret_sims <- (exp(cbind(EM_FTSE_ret_sims, EM_GSPC_ret_sims)) %*% port_weights) - 1

HF_VaR <- t(quantile(HF_port_ret_sims, c(0.01, 0.05)))
# EM_VaR <- t(quantile(EM_port_ret_sims, c(0.01, 0.05)))

# Set directory for VaR sims
save_dir <- paste0('./VaRBacktest/Sims')

# Write VaR sims to file
HF_VaR_filename <- paste0(save_dir, '/', task_id, '_HF_', family, '_VaR.csv')
# EM_VaR_filename <- paste0(save_dir, '/', task_id, '_EM_', family, '_VaR.csv')
write.table(HF_VaR, HF_VaR_filename, sep = ',', row.names = F, col.names = F)
# write.table(EM_VaR, EM_VaR_filename, sep = ',', row.names = F, col.names = F)

# Write copula GOF data to file
HF_copula_GOF <- t(c(HFResult$loglik, HFResult$aic, HFResult$bic))
# EM_copula_GOF <- t(c(EMResult$loglik, EMResult$aic, EMResult$bic))
HF_copula_GOF_filename <- paste0(save_dir, '/', task_id, '_HF_', family, '_GOF.csv')
# EM_copula_GOF_filename <- paste0(save_dir, '/', task_id, '_EM_', family, '_GOF.csv')
write.table(HF_copula_GOF, HF_copula_GOF_filename, sep = ',', row.names = F, col.names = F)
# write.table(EM_copula_GOF, EM_copula_GOF_filename, sep = ',', row.names = F, col.names = F)
      