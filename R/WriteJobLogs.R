WriteJobLogs <- function(job_options, log_dir) {
    names <- names(job_options)
    values <- numeric(length(job_options))
    for (i in 1:length(job_options)) {
      values[i] <- paste0(job_options[[i]], collapse = ', ')
    }
    
    job_options_mat <- cbind(names, values)
    write.csv(job_options_mat, paste0(log_dir, '/job', job_options$job_id, '_logs.csv'), row.names = F)
}