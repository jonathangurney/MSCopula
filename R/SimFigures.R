# Create plots of statistics from the data simulation and typical examples of
# simulated data

# Import packages
library(RColorBrewer)

job_id <- 1
family <- 'Clayton'
job_dir <- paste0('~/UCL MSc Statistics/STAT0034 MSc Project/Code/SimResults/', job_id)
data_dir <- paste0(job_dir)
fig_dir <- paste0('~/UCL MSc Statistics/STAT0034 MSc Project/Code/SimResults/Figures/', job_id)

colours <- brewer.pal(7, 'Pastel1')[c(1:5, 7)]

states_filename <- paste0(data_dir, '/', job_id, '_states.csv')
theta1_filename <- paste0(data_dir, '/', job_id, '_theta1.csv')
theta2_filename <- paste0(data_dir, '/', job_id, '_theta2.csv')
theta_sim_filename <- paste0(data_dir, '/', job_id, '_theta_sim.csv')

states <- as.matrix(read.csv(states_filename, header = F))
theta1 <- as.matrix(read.csv(theta1_filename, header = F))
theta2 <- as.matrix(read.csv(theta2_filename, header = F))
theta_sim <- as.matrix(read.csv(theta_sim_filename, header = F))

#
# Create histogram of dependence parameters in each state
#

lower <- min(c(theta1, theta2))
upper <- max(c(theta1, theta2))
range <- upper - lower
xlim <- c(lower - 0.05*range, upper + 0.05*range)
breaks <- seq(xlim[1], xlim[2], length.out = 50)

hist1 <- hist(theta1, breaks = breaks, plot = F)
hist2 <- hist(theta2, breaks = breaks, plot = F)

ylim <- c(0, 1.1*max(c(hist1$density, hist2$density)))

main = paste0('Histogram of copula dependence parameters for ', family, ' family')

cop_param_hist_filename <- paste0(fig_dir, '/', job_id, '_', family, '_dependence_param_hist.png')
png(filename = cop_param_hist_filename, width = 1200, height = 800, pointsize = 16)
par(cex = 1.25)
hist(theta1, breaks = breaks, xlim = xlim, ylim = ylim, freq = F, col = colours[1],
           main = main, xlab = expression(paste('Copula parameter, ', theta)))
hist(theta2, breaks = breaks, add = T, col = colours[2], freq = F)
legend('topright', c('State 1', 'State 2'), fill = colours[1:2], cex = 1.3,
       border = "black", box.lty = 'blank', inset = c(0.3, 0))
dev.off()

#
# Create plot of copula densities for typical state dynamics
#

# param_list <- c(6, 2, 5, 2, 0.9, 0.7)
# 
# clayton_cop1 <- copula::claytonCopula(param = param_list[1])
# clayton_cop2 <- copula::claytonCopula(param = param_list[2])
# 
# gumbel_cop1 <- copula::gumbelCopula(param = param_list[3])
# gumbel_cop2 <- copula::gumbelCopula(param = param_list[4])
# 
# gaussian_cop1 <- copula::normalCopula(param = param_list[5])
# gaussian_cop2 <- copula::normalCopula(param = param_list[6])
# 
# n_sims <- 3000
# 
# data <- list()
# 
# data[[1]] <- copula::rCopula(n_sims, clayton_cop1)
# data[[2]] <- copula::rCopula(n_sims, clayton_cop2)
# data[[3]] <- copula::rCopula(n_sims, gumbel_cop1)
# data[[4]] <- copula::rCopula(n_sims, gumbel_cop2)
# data[[5]] <- copula::rCopula(n_sims, gaussian_cop1)
# data[[6]] <- copula::rCopula(n_sims, gaussian_cop2)
# 
# fig_filename <- 'C:/Users/jonat/OneDrive/Documents/UCL MSc Statistics/STAT0034 MSc Project/Report figures/data_sim.png'
# png(fig_filename, width = 1200, height = 1200, pointsize = 16)
# par(mfrow = c(3, 2), cex.axis = 1.15, cex.lab = 1.5)
# for (i in 1:3) {
#   xlabel1 <- bquote(theta == .(param_list[2*i - 1]))
#   xlabel2 <- bquote(theta == .(param_list[2*i]))
#   plot(data[[2*i - 1]], pch = 19, xaxt = 'n', yaxt = 'n', col = colours[2*i - 1], bty = 'l', xlab = xlabel1, ylab = "")
#   axis(1, xaxp = c(0, 1, 1))
#   axis(2, yaxp = c(0, 1, 1))
#   plot(data[[2*i]], pch = 19, col = colours[2*i], xaxt = 'n', yaxt = 'n', bty = 'l', xlab = xlabel2, ylab = "")
#   axis(1, xaxp = c(0, 1, 1))
#   axis(2, yaxp = c(0, 1, 1))
# }
# mtext('Clayton copula', side = 3, line = -2, at = 0.5, padj = 1, cex = 1.5, outer = TRUE)
# mtext('Gumbel copula', side = 3, line = -32, at = 0.5, padj = 1, cex = 1.5, outer = TRUE)
# mtext('Gaussian copula', side = 3, line = -64, at = 0.5, padj = 1, cex = 1.5, outer = TRUE)
# dev.off()


