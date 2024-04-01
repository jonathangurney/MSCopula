# Transfer data files to another job with correct naming

setwd('~/Scratch/Code')

job_id_old <- 1
job_id_new <- 4

job_dir_old <- paste0('./SimResults/', job_id_old)
job_dir_new <- paste0('./SimResults/', job_id_new)

replications <- 2000

data_dir_old <- paste0(job_dir_old, '/Data')
data_dir_new <- paste0(job_dir_new, '/Data')

for (i in 1:replications) {
  data_filename_old <- paste0(data_dir_old, '/', job_id_old, '_', i, '_data.csv')
  data_filename_new <- paste0(data_dir_new, '/', job_id_new, '_', i, '_data.csv')
  
  states_filename_old <- paste0(data_dir_old, '/', job_id_old, '_', i, '_states.csv')
  states_filename_new <- paste0(data_dir_new, '/', job_id_new, '_', i, '_states.csv')
  
  theta_filename_old <- paste0(data_dir_old, '/', job_id_old, '_', i, '_theta.csv')
  theta_filename_new <- paste0(data_dir_new, '/', job_id_new, '_', i, '_theta.csv')
  
  theta_sim_filename_old <- paste0(data_dir_old, '/', job_id_old, '_', i, '_theta_sim.csv')
  theta_sim_filename_new <- paste0(data_dir_new, '/', job_id_new, '_', i, '_theta_sim.csv')
  
  file.copy(data_filename_old, data_filename_new)
  file.copy(states_filename_old, states_filename_new)
  file.copy(theta_filename_old, theta_filename_new)
  file.copy(theta_sim_filename_old, theta_sim_filename_new)
}
