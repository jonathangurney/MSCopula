# Collate VaR backtest results

CollateVaRBacktest <- function() {
  n_dates <- 1470
  families <- c('clayton', 'gumbel', 'gaussian')
  
  data_dir <- '~/Scratch/Code/VaRBacktest/Sims'
  
  HF_VaR <- matrix(nrow = n_dates, ncol = 6)
  EM_VaR <- matrix(nrow = n_dates, ncol = 6)
  HF_GOF <- matrix(nrow = n_dates, ncol = 9)
  EM_GOF <- matrix(nrow = n_dates, ncol = 9)
  
  for (i in 1:n_dates) {
    if (i %% 50 == 0) {
      print(paste0(i, '/', n_dates))
    }
    for (family in families) {
      j <- which(families == family)
      
      HF_VaR_data_filename <- paste0(data_dir, '/', i, '_HF_', family, '_VaR.csv')
      EM_VaR_data_filename <- paste0(data_dir, '/', i, '_EM_', family, '_VaR.csv')
      HF_GOF_data_filename <- paste0(data_dir, '/', i, '_HF_', family, '_GOF.csv')
      EM_GOF_data_filename <- paste0(data_dir, '/', i, '_EM_', family, '_GOF.csv')
      
      if (!file.exists(HF_VaR_data_filename)) {
        print(paste0(i, ': HF data missing for ', family, ' copula'))
        next
      }
      if (!file.exists(EM_VaR_data_filename)) {
        print(paste0(i, ': EM data missing for ', family, ' copula'))
      }
      
      HF_VaR_data <- as.matrix(read.csv(HF_VaR_data_filename, header = F))
      EM_VaR_data <- as.matrix(read.csv(EM_VaR_data_filename, header = F))
      HF_GOF_data <- as.matrix(read.csv(HF_GOF_data_filename, header = F))
      EM_GOF_data <- as.matrix(read.csv(EM_GOF_data_filename, header = F))
      
      HF_VaR[i, (2*j - 1):(2*j)] <- HF_VaR_data
      EM_VaR[i, (2*j - 1):(2*j)] <- EM_VaR_data
      HF_GOF[i, (3*j - 2):(3*j)] <- HF_GOF_data
      EM_GOF[i, (3*j - 2):(3*j)] <- EM_GOF_data
    }
  }
  
  colnames(HF_VaR) <- as.character(t(outer(families, c(0.99, 0.95), paste)))
  colnames(EM_VaR) <- as.character(t(outer(families, c(0.99, 0.95), paste)))
  colnames(HF_GOF) <- as.character(t(outer(families, c('LL', 'AIC', 'BIC'), paste)))
  colnames(EM_GOF) <- as.character(t(outer(families, c('LL', 'AIC', 'BIC'), paste)))
  
  results_dir <- '~/Scratch/Code/VaRBacktest/Results'
  
  HF_VaR_filename <- paste0(results_dir, '/HF_VaR_VaR.csv')
  EM_VaR_filename <- paste0(results_dir, '/EM_VaR_VaR.csv')
  HF_GOF_filename <- paste0(results_dir, '/HF_VaR_GOF.csv')
  EM_GOF_filename <- paste0(results_dir, '/EM_VaR_GOF.csv')
  
  write.csv(HF_VaR, HF_VaR_filename, row.names = F)
  write.csv(EM_VaR, EM_VaR_filename, row.names = F)
  write.csv(HF_GOF, HF_GOF_filename, row.names = F)
  write.csv(EM_GOF, EM_GOF_filename, row.names = F)
}

CollateVaRBacktest()