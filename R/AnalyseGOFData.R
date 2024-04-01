# Set working directory
setwd('~/UCL MSc Statistics/STAT0034 MSc Project/Code')

# Import packages
library(RColorBrewer)

# Set path for simulation results
sim_results_path <- './SimResults'

# Create vector declaring copula families
job_id <- 6
family <- 'gaussian'
family_nice <- 'Gaussian'

families <- c('clayton', 'gumbel', 'gaussian')
families_nice <- c('Clayton', 'Gumbel', 'Gaussian')
n_families <- length(families)

observations = 1000
replications <- 2000

data_dir <- paste0(sim_results_path, '/', job_id, '/DataFinal') 

# Initialize lists for GOF data
HF_GOF <- list()
EM_GOF <- list()
for (f in families) {
  HF_GOF[[f]] <- read.csv(paste0(data_dir, '/', 'HF_', f, '_GOF.csv'), header = F)
  EM_GOF[[f]] <- read.csv(paste0(data_dir, '/', 'EM_', f, '_GOF.csv'), header = F)
}
# Create plots

# Set colour palette
colours <- brewer.pal(name = 'Pastel1', n = 5)

# First check if correct figure directory exists
fig_dir <- paste0(sim_results_path, '/Figures/', job_id)
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}

gof_types <- c("Log-likelihood", "AIC", "BIC")

# First plot will compare Log-Likelihood, AIC and BIC distributions between
# estimation methods
gof_fig_filename <- paste0(fig_dir, '/gof_plot_1.png')
png(filename = gof_fig_filename, width = 1200, height = 800, pointsize = 16)
par(mfrow = c(1, 3), cex = 1.15)
for (i in 1:3) {
  boxplot(list(HF_GOF[[family]][ , i], EM_GOF[[family]][ , i]),
          xlab = gof_types[i],
          names = c("IFM", "EM"),
          col = colours[1:2],
          pch = 19, cex = 0.8, lwd = 1.2, frame.plot = F, whisklwd = 2)
}
mtext(paste0("GOF distributions for ", family_nice, " copula with n = ", observations),
      side = 3, line = -2.5, outer = TRUE, cex = 1.8)
dev.off()

# Second plot will compare BIC values for all copula families fitted to each data set

gof_fig_2_filename <- paste0(fig_dir, '/gof_plot_2.png')
png(filename = gof_fig_2_filename, width = 1200, height = 800, pointsize = 16)
par(mfrow = c(1, 3), cex = 1.15)
for (i in 1:3) {
  par(cex.lab = 1.15, cex.axis = 1.15)
  boxplot(list(HF_GOF[[families[i]]][ , 3], EM_GOF[[families[i]]][ , 3]),
          xlab = families_nice[i],
          ylab = 'BIC',
          names = c('IFM', 'EM'),
          col = colours[1:2],
          pch = 19, cex = 0.8, lwd = 1.2, frame.plot = F, whisklwd = 2)
  box(bty = 'l')
}
dev.off()

# Create table for how often each model is fitted to the data based on BIC criterion

model_selected <- matrix(nrow = replications, ncol = 2)
for (i in 1:replications) {
  HF_bic <- numeric(n_families)
  EM_bic <- numeric(n_families)
  for (j in 1:n_families) {
    HF_bic[j] <- HF_GOF[[families[j]]][i, 3]
    EM_bic[j] <- EM_GOF[[families[j]]][i, 3]
  }

  if (!any(is.na(HF_bic))) {
    model_selected[i, 1] <- which(HF_bic == min(HF_bic))
  }
  if (!any(is.na(EM_bic))) {
    model_selected[i, 2] <- which(EM_bic == min(EM_bic))
  }
}
HF_GOF_table <- table(model_selected[ , 1])
EM_GOF_table <- table(model_selected[ , 2])

print(100*prop.table(HF_GOF_table))
print(100*prop.table(EM_GOF_table))
