# Write a function to collate the results of the n_sims into one file per
# important category

CollateResults <- function(jobs, families) {
  n_obs <- 1000
  n_sims <- 2000
  
  for (job_idx in jobs) {
    job_dir <- paste0('./SimResults/', job_idx)
    data_dir <- paste0(job_dir, '/Data')
    
    for (family in families) {
      print(family)
      
      # Initialize the matrices for storing the data
      HF_params <- matrix(nrow = n_sims, ncol = 6)
      HF_GOF <- matrix(nrow = n_sims, ncol = 3)
      HF_errors <- matrix(nrow = n_sims, ncol = 3)
      # EM_params <- matrix(nrow = n_sims, ncol = 6)
      # EM_GOF <- matrix(nrow = n_sims, ncol = 3)
      # EM_errors <- matrix(nrow = n_sims, ncol = 3)
      
      for (i in 1:n_sims) {
        if ((i %% 50) == 0) {
          print(paste0(i, '/', n_sims))
        }
        
        HF_params_filename <- paste0(data_dir, '/', job_idx, '_', i, '_', family, '_HF_parameters_.csv')
        HF_GOF_filename <- paste0(data_dir, '/', job_idx, '_', i, '_', family, '_HF_GOF_.csv')
        HF_errors_filename <- paste0(data_dir, '/', job_idx, '_', i, '_', family, '_HF_errors_.csv')
        # EM_params_filename <- paste0(data_dir, '/', job_idx, '_', i, '_', family, '_EM_parameters_.csv')
        # EM_GOF_filename <- paste0(data_dir, '/', job_idx, '_', i, '_', family, '_EM_GOF_.csv')
        # EM_errors_filename <- paste0(data_dir, '/', job_idx, '_', i, '_', family, '_EM_errors_.csv')
        
        # Check if the necessary files exist
        HF_params_exists <- file.exists(HF_params_filename)
        HF_GOF_exists <- file.exists(HF_GOF_filename)
        HF_errors_exists <- file.exists(HF_errors_filename)
        # EM_params_exists <- file.exists(EM_params_filename)
        # EM_GOF_exists <- file.exists(EM_GOF_filename)
        # EM_errors_exists <- file.exists(EM_errors_filename)
        
        if (!all(HF_params_exists, HF_GOF_exists, HF_errors_exists)) {
                 # EM_params_exists, EM_GOF_exists, EM_errors_exists)) {
          print(paste0(job_idx, ': Files do not exist for ', family, ' copula simulation ', i))
          next
        }
        
        HF_params_data <- as.matrix(read.csv(HF_params_filename, header = F))
        HF_GOF_data <- as.matrix(read.csv(HF_GOF_filename, header = F))
        HF_errors_data <- as.matrix(read.csv(HF_errors_filename, header = F))
        # EM_params_data <- as.matrix(read.csv(EM_params_filename, header = F))
        # EM_GOF_data <- as.matrix(read.csv(EM_GOF_filename, header = F))
        # EM_errors_data <- as.matrix(read.csv(EM_errors_filename, header = F))
        
        HF_params[i, ] <- HF_params_data
        
        # Make minor adjustment for when AIC, BIC hasn't been calculated for some
        # models
        if (ncol(HF_GOF_data) == 1) {
          aic <- 2 * 6 - 2 * HF_GOF_data[1]
          bic <- 6 * log(2 * n_obs) - 2 * HF_GOF_data[1]
          HF_GOF_data <- c(HF_GOF_data[1], aic, bic)
        }
        
        HF_GOF[i, ] <- HF_GOF_data
        HF_errors[i, ] <- HF_errors_data
        # EM_params[i, ] <- EM_params_data
        # EM_GOF[i, ] <- EM_GOF_data
        # EM_errors[i, ] <- EM_errors_data
      }
      
      save_data_dir <- paste0(job_dir, '/DataFinal')
      if (!dir.exists(save_data_dir)) {
        dir.create(save_data_dir)
      }
      
      HF_params_save_file <- paste0(save_data_dir, '/HF_', family, '_params.csv')
      HF_GOF_save_file <- paste0(save_data_dir, '/HF_', family, '_GOF.csv')
      HF_errors_save_file <- paste0(save_data_dir, '/HF_', family, '_errors.csv')
      # EM_params_save_file <- paste0(save_data_dir, '/EM_', family, '_params.csv')
      # EM_GOF_save_file <- paste0(save_data_dir, '/EM_', family, '_GOF.csv')
      # EM_errors_save_file <- paste0(save_data_dir, '/EM_', family, '_errors.csv')
      
      write.table(HF_params, HF_params_save_file, row.names = F, col.names = F, sep = ',')
      write.table(HF_GOF, HF_GOF_save_file, row.names = F, col.names = F, sep = ',')
      write.table(HF_errors, HF_errors_save_file, row.names = F, col.names = F, sep = ',')
      # write.table(EM_params, EM_params_save_file, row.names = F, col.names = F, sep = ',')
      # write.table(EM_GOF, EM_GOF_save_file, row.names = F, col.names = F, sep = ',')
      # write.table(EM_errors, EM_errors_save_file, row.names = F, col.names = F, sep = ',')
    }
  }
}

setwd('~/Scratch/Code')
job_idx <- c(4, 5, 6)
families <- c('clayton', 'gaussian', 'gumbel')

CollateResults(jobs = job_idx, families = families)
