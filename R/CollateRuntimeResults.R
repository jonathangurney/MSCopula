CollateRuntimeResults <- function(jobs, families, include_EM) {
  n_obs <- 1000
  n_sims <- 2000
  
  for (job_id in jobs) {
    job_index <- which(jobs == job_id)
    
    job_dir <- paste0('./SimResults/', job_id)
    data_dir <- paste0(job_dir, '/Data')
    
    for (family in families) {
      print(family)
      
      print(job_index)
      print(include_EM[job_index])
      
      # Initialize matrix for storing data
      if (include_EM[job_index]) {
        runtimes <- matrix(nrow = n_sims, ncol = 2)
      } else {
        runtimes <- matrix(nrow = n_sims, ncol = 1)
      }
        
      for (i in 1:n_sims) {
        if ((i %% 50) == 0) {
          print(paste0(i, '/', n_sims))
        }
        
        runtimes_data_filename <- paste0(data_dir, '/', job_id, '_', i, '_', family, '_runtimes_.csv')
        
        if (!file.exists(runtimes_data_filename)) {
          print(paste0(job_id, ': Runtimes file does not exist for simulation ', i))
          next
        }
        
        runtimes_data <- as.matrix(read.csv(runtimes_data_filename, header = F))
        runtimes[i , ] <- runtimes_data
      }
     
      save_data_dir <- paste0(job_dir, '/DataFinal')
      if (!dir.exists(save_data_dir)) {
        dir.create(save_data_dir)
      }
      
      runtimes_save_filename <- paste0(save_data_dir, '/', job_id, '_', family, '_runtimes.csv')
                                       
      if (include_EM[job_index]) {
        colnames(runtimes) <- c('HF_runtime', 'EM_runtime')
        write.csv(runtimes, runtimes_save_filename, row.names = F)
      } else {
        colnames(runtimes) <- c('HF_runtime')
        write.csv(runtimes, runtimes_save_filename, row.names = F)
      }
      
    }
  }
}
      
setwd('~/Scratch/Code')
job_idx <- c(1, 2, 3, 4, 5, 6)
include_EM <- c(T, T, T, F, F, F)
families <- c('clayton', 'gaussian', 'gumbel')

CollateRuntimeResults(jobs = job_idx, families = families, include_EM = include_EM)