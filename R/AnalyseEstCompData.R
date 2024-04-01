# Set working directory
setwd('~/UCL MSc Statistics/STAT0034 MSc Project/Code')

# Import packages
library(RColorBrewer)

# Set path for simulation results
sim_results_path <- './SimResults'

# Create vector declaring copula families
job_idx <- 4
family <- 'clayton'
family_nice <- 'Clayton'
print(family_nice)

observations = 1000

data_dir <- paste0(sim_results_path, '/', job_idx, '/DataFinal')

# Read csv files
HFErrors <- read.csv(paste0(data_dir, '/', 'HF_', family, '_errors.csv'), header = F)
EMErrors <- read.csv(paste0(data_dir, '/', 'EM_', family, '_errors.csv'), header = F)
HFPars <- read.csv(paste0(data_dir, '/', 'HF_', family, '_params.csv'), header = F)
EMPars <- read.csv(paste0(data_dir, '/', 'EM_', family, '_params.csv'), header = F)

# Create plots

# Set colour palette
colours <- brewer.pal(5, 'Pastel1')

# First check if correct figure directory exists
fig_dir <- paste0(sim_results_path, '/Figures/', job_idx)
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}

# error_types <- c("Mean Error", "RMSE", "MAE")
# arma_parameters <- c(expression(omega[1]), expression(omega[2]), expression(paste(alpha)), expression(paste(beta)))
# 
# arma_pars_list <- list(gaussian = c(2.5, 1, 0.3, 0.4),
#                   gumbel = c(1.7, 0.7, 0.1, -1),
#                   clayton = c(2, 1.4, 0.1, -1),
#                   frank = c(3, 1.5, 0.3, -1.1))
# 
# # First plot will compare errors in theta across estimation methods
# err_fig_filename <- paste0(fig_dir, '/error_plot_', observations, '.png')
# png(filename = err_fig_filename, width = 1200, height = 800, pointsize = 16)
# par(mfrow = c(1, 3), cex = 1.15)
# for (i in 1:3) {
#   boxplot(list(HFErrors[ , i], EMErrors[ , i]),
#           xlab = error_types[i],
#           names = c("IFM", "EM"),
#           col = colours[1:2],
#           pch = 19, cex = 0.8, lwd = 1.2, frame.plot = F, whisklwd = 2)
#   box(bty = 'l')
# }
# mtext(paste0("Error distributions for ", family_nice, " copula with n = ", observations),
#       side = 3, line = -2.5, outer = TRUE, cex = 1.8)
# dev.off()

# Create a table with summary statistics for the dependence parameter error distributions

HF_theta_errors <- list()
EM_theta_errors <- list()

for (i in 1:3) {
  HF_theta_errors[[i]] <- summary(HFErrors[ , i])
  EM_theta_errors[[i]] <- summary(EMErrors[ , i])
}
stats <- c('Min.', 'Median', 'Mean', 'Max.')
stats_nice <- c('Minimum', 'Median', 'Mean', 'Maximum')

cat(paste0(' & \\multicolumn{2}{c|}{Mean Error} & \\multicolumn{2}{c|}{RMSE} & \\multicolumn{2}{c}{MAE} \\\\\n'))
cat(paste0(' & ', paste0(rep(c('\\textit{IFM}', '\\textit{EM}'), 3), collapse = ' & '), '\\\\\n'))
cat('\\hline\n')
for (i in 1:4) {
  cat(paste0(stats_nice[i]))
  for (j in 1:3) {
    cat(paste0(' & ', formatC(HF_theta_errors[[j]][stats[i]], digits = 3, format = 'f'), ' & '))
    cat(paste0(formatC(EM_theta_errors[[j]][stats[i]], digits = 3, format = 'f')))
  }
  cat(' \\\\\n')
}
cat('Standard deviation')
for (j in 1:3) {
  cat(paste0(' & ', formatC(sd(HFErrors[ , j], na.rm = T), digits = 3, format = 'f'), ' & '))
  cat(paste0(formatC(sd(EMErrors[ , j], na.rm = T), digits = 3, format = 'f')))
}
cat(' \\\\\n')

# Now look to plot parameter estimates produced by each estimation procedure

# We can break this down into two plots - one for the parameters for the ARMA
# process and one for the transition matrix probabilities

# Need to modify the data so that we consider the regimes correctly. Then we can
# consider transition probabilities for the high and low dependence regimes
# HFParsSorted <- HFPars
# EMParsSorted <- EMPars
# print(sum(HFPars[ , 1] < HFPars[ , 2]))
# print(sum(EMPars[ , 1] < HFPars[ , 2]))
# for (i in 1:nrow(HFPars)) {
#   if (any(is.na(HFPars[i, ]))) {
#     next
#   }
# 
#   if (any(is.na(EMPars[i, ]))) {
#     next
#   }
# 
#   if (HFPars[i, 1] < HFPars[i, 2]) {
#     HFParsSorted[i, 1] <- HFPars[i, 2]
#     HFParsSorted[i, 2] <- HFPars[i, 1]
#     HFParsSorted[i, 3:4] <- HFPars[i, 3:4]
#     HFParsSorted[i, 5] <- -HFPars[i, 6]
#     HFParsSorted[i, 6] <- -HFPars[i, 5]
#   }
# 
#   if (EMPars[i, 1] < EMPars[i, 2]) {
#     EMParsSorted[i, 1] <- EMPars[i, 2]
#     EMParsSorted[i, 2] <- EMPars[i, 1]
#     EMParsSorted[i, 3:4] <- EMPars[i, 3:4]
#     EMParsSorted[i, 5] <- -EMPars[i, 6]
#     EMParsSorted[i, 6] <- -EMPars[i, 5]
#   }
# }
# # 
# # arma_pars_fig_filename <- paste0(fig_dir, '/arma_pars_plot_sorted_', observations, '.png')
# # png(filename = arma_pars_fig_filename, width = 1200, height = 800, pointsize = 16)
# # par(mfrow = c(1, 4), cex = 1.15, cex.lab = 1.5)
# # for (i in 1:4) {
# #   boxplot(list(HFParsSorted[ , i], EMParsSorted[ , i]),
# #           xlab = arma_parameters[i],
# #           names = c("IFM", "EM"),
# #           col = colours[1:2],
# #           pch = 19, cex = 0.8, lwd = 1.2, frame.plot = F, whisklwd = 2)
# #   abline(h = arma_pars_list[[family]][i], col = 'red', lwd = 2)
# #   box(bty = 'l')
# # }
# # mtext(paste0("ARMA parameter estimates for ", family_nice, " copula with n = ", observations, ' observations'), side = 3, line = -2.5, outer = TRUE, cex = 1.8)
# # dev.off()
# # 
# # # Create plot of transition probabilities
# # 
# # trans_probs <- c(0.9, 0.8)
# # 
# # # Define the logistic function for transforming the parameters to probabilities
# # logistic <- function(x) {
# #   return(1/(1.00001 + exp(-x)))
# # }
# # 
# # trans_pars_fig_filename <- paste0(fig_dir, '/trans_pars_plot_sorted_', observations, '.png')
# # png(filename = trans_pars_fig_filename, width = 1200, height = 800, pointsize = 16)
# # par(mfrow = c(1, 2), cex = 1.15, cex.lab = 1.2)
# # boxplot(list(logistic(HFParsSorted[ , 5]), logistic(EMParsSorted[ , 5])),
# #         xlab = expression(P[11]),
# #         names = c("IFM", "EM"),
# #         col = colours[1:2],
# #         pch = 19, cex = 0.8, lwd = 1.2, frame.plot = F, whisklwd = 2, ylim = c(0, 1))
# # abline(h = trans_probs[1], col = 'red', lwd = 2)
# # box(bty = 'l')
# # boxplot(list((1 - logistic(HFParsSorted[ , 6])), (1 - logistic(EMParsSorted[ , 6]))),
# #         xlab = expression(P[22]),
# #         names = c("IFM", "EM"),
# #         col = colours[1:2],
# #         pch = 19, cex = 0.8, lwd = 1.2, frame.plot = F, whisklwd = 2, ylim = c(0, 1))
# # abline(h = trans_probs[2], col = 'red', lwd = 2)
# # box(bty = 'l')
# # mtext(paste0("Transition probability parameter estimates for ", family_nice, " copula with n = ", observations), side = 3, line = -2.5, outer = TRUE, cex = 1.8)
# # dev.off()
# # 
# # # Compute summary statistics for the ARMA parameter distributions
# # 
# # stat_names = c('Minimum', 'Median', 'Mean', 'Maximum', 'Standard Deviation')
# # 
# # HF_param_stats <- vector('list', 4)
# # EM_param_stats <- vector('list', 4)
# # for (i in 1:4) {
# #   HF_param_stats[[i]] <- c(summary(HFParsSorted[ , i] - arma_pars_list[[family]][i])[c('Min.', 'Median', 'Mean', 'Max.')], SD = sd(HFParsSorted[ , i], na.rm = T))
# #   EM_param_stats[[i]] <- c(summary(EMParsSorted[ , i] - arma_pars_list[[family]][i])[c('Min.', 'Median', 'Mean', 'Max.')], SD = sd(EMParsSorted[ , i], na.rm = T))
# # }
# # 
# # cat(' & \\multicolumn{2}{|c|}{$\\omega_1$} & \\multicolumn{2}{|c|}{$\\omega_2$} & \\multicolumn{2}{|c|}{$\\alpha$} & \\multicolumn{2}{|c|}{$\\beta$} \\\\\n')
# # cat(paste0(' & ', paste0(rep(c('\\textit{IFM}', '\\textit{EM}'), 4), collapse = ' & '), '\\\\\n'))
# # cat('\\hline\n')
# # for (i in 1:5) {
# #   cat(paste0(stat_names[i], ' & '))
# #   for (j in 1:4) {
# #     cat(paste0(formatC(HF_param_stats[[j]][i], digits = 4, format = 'f'), ' & ',
# #                formatC(EM_param_stats[[j]][i], digits = 4, format = 'f'), ' & '))
# #   }
# #   cat('\\\\\n')
# # }
# 
# # Compute summary statistics for transition probabilities
# 
# logistic <- function(x) {
#   return(1/(1 + exp(-x)))
# }
# 
# HF_p11_stats <- summary(logistic(HFParsSorted[ , 5]) - 0.9)
# HF_p22_stats <- summary(logistic(-HFParsSorted[ , 6]) - 0.8)
# EM_p11_stats <- summary(logistic(EMParsSorted[ , 5]) - 0.9)
# EM_p22_stats <- summary(logistic(-EMParsSorted[ , 6]) - 0.8)
# stats <- c('Min.', 'Median', 'Mean', 'Max.')
# stats_nice <- c('Minimum', 'Median', 'Mean', 'Maximum')
# 
# cat(paste0(' & \\multicolumn{2}{c|}{$P_{11}$} & \\multicolumn{2}{c|}{$P_{22}$} \\\\\n'))
# cat(paste0(' & ', paste0(rep(c('\\textit{IFM}', '\\textit{EM}'), 2), collapse = ' & '), '\\\\\n'))
# cat('\\hline\n')
# for (i in 1:4) {
#   cat(paste0(stats_nice[i], ' & '))
#   cat(paste0(formatC(HF_p11_stats[stats[i]], digits = 3, format = 'f'), ' & '))
#   cat(paste0(formatC(EM_p11_stats[stats[i]], digits = 3, format = 'f'), ' & '))
#   cat(paste0(formatC(HF_p22_stats[stats[i]], digits = 3, format = 'f'), ' & '))
#   cat(paste0(formatC(EM_p22_stats[stats[i]], digits = 3, format = 'f'), ' \\\\\n'))
# }
# cat('Standard deviation & ')
# cat(paste0(formatC(sd(logistic(HFParsSorted[ , 5]), na.rm = T), digits = 3, format = 'f'), ' & '))
# cat(paste0(formatC(sd(logistic(EMParsSorted[ , 5]), na.rm = T), digits = 3, format = 'f'), ' & '))
# cat(paste0(formatC(sd(logistic(-HFParsSorted[ , 6]), na.rm = T), digits = 3, format = 'f'), ' & '))
# cat(paste0(formatC(sd(logistic(-EMParsSorted[ , 6]), na.rm = T), digits = 3, format = 'f'), ' \\\\\n'))
# 
# 
# 
# 
# 
