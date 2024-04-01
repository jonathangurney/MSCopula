# Aggregation of state statistics, theta values etc

CollateSimData <- function(job_ids) {
  for (job_id in job_ids) {
    job_dir <- paste0('~/Scratch/Code/SimResults/', job_id)
    data_dir <- paste0(job_dir, '/Data')
    data_final_dir <- paste0(job_dir, '/DataFinal')
    
    observations <- 1000
    replications <- 2000
    
    # Initialize dataframes for storing collated data
    states_data <- matrix(nrow = observations, ncol = replications)
    theta1_data <- matrix(nrow = observations, ncol = replications)
    theta2_data <- matrix(nrow = observations, ncol = replications)
    theta_sim_data <- matrix(nrow = observations, ncol = replications)
    
    for (i in 1:replications) {
      if ((i %% 100) == 0) {
        print(paste0(i, '/', replications))
      }

      states_filename <- paste0(data_dir, '/', job_id, '_', i, '_states.csv')
      theta_filename <- paste0(data_dir, '/', job_id, '_', i, '_theta.csv')
      theta_sim_filename <- paste0(data_dir, '/', job_id, '_', i, '_theta_sim.csv')
      
      state <- as.matrix(read.csv(states_filename, header = F))
      theta <- as.matrix(read.csv(theta_filename, header = F))
      theta_sim <- as.matrix(read.csv(theta_sim_filename, header = F))
      
      states_data[ , i] <- state
      theta1_data[ , i] <- theta[ , 1]
      theta2_data[ , i] <- theta[ , 2]
      theta_sim_data[ , i] <- theta_sim
    }
    
    states_data_filename <- paste0(data_final_dir, '/', job_id, '_states.csv')
    theta1_data_filename <- paste0(data_final_dir, '/', job_id, '_theta1.csv')
    theta2_data_filename <- paste0(data_final_dir, '/', job_id, '_theta2.csv')
    theta_sim_data_filename <- paste0(data_final_dir, '/', job_id, '_theta_sim.csv')
    
    write.table(states_data, states_data_filename, row.names = F, col.names = F, sep = ',')
    write.table(theta1_data, theta1_data_filename, row.names = F, col.names = F, sep = ',')
    write.table(theta2_data, theta2_data_filename, row.names = F, col.names = F, sep = ',')
    write.table(theta_sim_data, theta_sim_data_filename, row.names = F, col.names = F, sep = ',')
  }
}

jobs <- c(1, 2, 3)

CollateSimData(jobs)

