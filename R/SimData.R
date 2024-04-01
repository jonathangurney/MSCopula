SimData <- function(n, pars, regimes = 2, trans_mat, family,
                    burn_in = 1000, save = FALSE, dir = NA, ...) {
  params <- list(...)
  if (family == 't') {
    df <- params[['df']]
  }
  
  states_sim <- CreateMarkovChain((n + burn_in),
                                  regimes = regimes, trans_mat = trans_mat)
  
  # Establish the number of time varying parameters
  if (family == 'sjc') {
    p <- 2
  } else {
    p <- 1
  }
  
  # Determine the form of the psi terms
  if (family == 'gaussian') {
    psi_func <- function(u1, u2) {
      return(mean(qnorm(u1) * qnorm(u2)))
    }
  } else if (family == 't') {
    psi_func <- function(u1, u2, df = 4) {
      return(mean(qt(u1, df = df) * qt(u2, df = df)))
    }
  } else {
    psi_func <- function(u1, u2) {
      return(mean(abs(u1 - u2)))
    }
  }
  
  # Initialize
  theta <- vector("list", p)
  theta_sim <- matrix(nrow = (n + burn_in), ncol = p)
  
  switch(family,
         "gaussian" = {
           transform_func <- function(x) {
             return(1.9998/(exp(-x) + 1) - 0.9999)
           }
           
           theta_init <- 0.5
         },
         
         "gumbel" = {
           transform_func <- function(x) {
             return(1 + x^2)
           }
           
           theta_init <- 2
         },
         
         "clayton" = {
           transform_func <- function(x) {
             return(x^2)
           }
           theta_init <- 2
         },
         
         "frank" = {
           transform_func <- function(x) {
             return(x)
           }
           
           theta_init <- 2
         },
         
         "t" = {
           transform_func <- function(x) {
             return(1.9998/(exp(-x) + 1) - 0.9999)
           }
           
           theta_init <- 0.5
         })
  
  for (j in 1:p) {
    theta[[j]] <- matrix(nrow = (n + burn_in), ncol = regimes)
    theta[[j]][1, ] <- rep(theta_init, regimes)
  }
  
  theta_sim[1, ] <- rep(theta_init, p)
  
  U <- matrix(nrow = (n + burn_in), ncol = 2)
  
  if (family == 't') {
    U[1, ] <- cop_sim(1, theta_init, family, df = df[states_sim[1]])
  } else {
    U[1, ] <- cop_sim(1, theta_init, family)
  }
  
  for (i in 2:(n + burn_in)) {
    if (i <= 10) {
      psi <- psi_func(U[1:(i - 1), 1], U[1:(i - 1), 2])
    } else {
      psi <- psi_func(U[(i - 10):(i - 1), 1], U[(i - 10):(i - 1), 2])
    }
    
    for (j in 1:p) {
      k <- (j - 1) * (regimes + 2)
      lp <- pars[(k + 1):(k + regimes)] + pars[k + regimes + 1] *
        theta[[p]][(i - 1), ] + pars[k + regimes + 2] * psi
      theta[[p]][i, ] <- transform_func(lp)
      theta_sim[i, p] <- theta[[p]][i, states_sim[i]]
    }
    
    if (family == 't') {
      U[i, ] <- cop_sim(1, theta_sim[i, ], family, df = df[states_sim[i]])
    } else {
      U[i, ] <- cop_sim(1, theta_sim[i, ], family)
    }
  }
  
  for (j in 1:p) {
    theta[[j]] <- theta[[j]][(burn_in + 1):(burn_in + n), ]
  }
  
  output <- list(U = U[(burn_in + 1):(burn_in + n), ],
                 states = states_sim[(burn_in + 1):(burn_in + n)],
                 theta = theta,
                 theta_sim = theta_sim[(burn_in + 1):(burn_in + n), ])
  
  return(output)
}

cop_sim <- function(n, param, family, ...) {
  param <- as.numeric(param)
  add_params <- list(...)
  df <- add_params$df
  switch(family,
         "gaussian" = {
           cop <- copula::normalCopula(param = param)
         },
         
         "gumbel" = {
           cop <- copula::gumbelCopula(param = param)
         },
         
         "clayton" = {
           cop <- copula::claytonCopula(param = param)
         },
         
         "frank" = {
           cop <- copula::frankCopula(param = param)
         },
         
         "t" = {
           cop <- copula::tCopula(param = param, df = df)
         }
         )
  data <- copula::rCopula(n, cop)
  return(data)
}

CreateMarkovChain <- function(n, regimes, trans_mat) {
  states <- numeric(n)
  
  # Create initial state
  init_probs <- rep(1/regimes, regimes)
  states[1] <- sample(1:regimes, 1, prob = init_probs)
  
  for (i in 2:n) {
    state0 <- rep(0, regimes)
    state0[states[i - 1]] <- 1
    
    state1_probs <- trans_mat %*% state0
    states[i] <- sample(1:regimes, 1, prob = state1_probs)
  }
  
  return(states)
}

