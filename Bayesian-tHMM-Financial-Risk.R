library(nimble)
library(coda)
library(MASS)
library(mvtnorm)
library(quantmod)
library(ggplot2)
library(tidyr)
library(dplyr)
library(corrplot)
library(gridExtra)



get_stock_data <- function(symbols, start_date, end_date) {
  stock_data <- list()
  
  for (symbol in symbols) {
    tryCatch({
      data <- getSymbols(symbol, src = "yahoo", 
                         from = start_date, 
                         to = end_date, 
                         auto.assign = FALSE)
      stock_data[[symbol]] <- Cl(data)  
    })
  }
  return(stock_data)
}

symbols <- c("NVDA", "MSFT", "META", "AAPL", "GOOG", "^GSPC", "^NDX")
start_date <- "2020-01-01"
end_date <- "2025-12-01"

stock_prices <- get_stock_data(symbols, start_date, end_date)

log_returns_list <- lapply(stock_prices, function(p) {
  100 * diff(log(p))  
})

log_returns <- do.call(merge, log_returns_list)

colnames(log_returns) <- c("NVIDIA", "MICROSOFT", "META", "APPLE", "GOOGLE", "SEP500", "NASDAQ")

log_returns <- na.omit(log_returns)

dates_clean <- index(log_returns)         
Y <- coredata(log_returns)                


print_descriptive_stats <- function(data) {
  stats <- data.frame(
    Mean = colMeans(data),
    Std.Dev = apply(data, 2, sd),
    Skewness = apply(data, 2, function(x) moments::skewness(x)),
    Kurtosis = apply(data, 2, function(x) moments::kurtosis(x)),
    Min = apply(data, 2, min),
    Median = apply(data, 2, median),
    Max = apply(data, 2, max)
  )
  return(stats)
}

print(round(print_descriptive_stats(Y), 3))

correlation <- cor(Y)
corrplot(correlation, method = 'number')


original_symbols <- names(stock_prices)   
display_names    <- colnames(log_returns)  

for (i in seq_along(original_symbols)) {
  
  sym <- original_symbols[i]
  name <- display_names[i]
  
  df_price <- data.frame(
    Date = index(stock_prices[[sym]]),
    Value = as.numeric(stock_prices[[sym]]),
    Type = "Price ($)" 
  )
  
  df_return <- data.frame(
    Date = index(log_returns[, i]),
    Value = as.numeric(log_returns[, i]),
    Type = "Log Returns (%)" 
  )
  
  plot_data <- rbind(df_price, df_return)
  
  plot_data$Type <- factor(plot_data$Type, levels = c("Price ($)", "Log Returns (%)"))
  
  p <- ggplot(plot_data, aes(x = Date, y = Value)) +
    
    geom_line(color = "#2c3e50", size = 0.3) +
    
    facet_grid(Type ~ ., scales = "free_y", switch = "y") +
    
    theme_gray() +
    
    labs(title = paste(name," Analysis"), x = NULL, y = NULL) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      strip.background = element_rect(fill = "lightgrey"), 
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),   
      panel.grid.major = element_line(color = "white"),    
      panel.grid.minor = element_line(color = "white")
    )
  
  print(p)
  
  Sys.sleep(1)
}

##Nimble implementation

tHMM_code <- nimbleCode({
  
  # Priors for initial state probabilities
  pi[1:K] ~ ddirch(alpha[1:K])
  
  # Priors for transition probabilities
  for(u in 1:K) {
    Pi[u, 1:K] ~ ddirch(alpha[1:K])
  }
  
  # Priors for state-specific means
  for(u in 1:K) {
    mu[u, 1:D] ~ dmnorm(m_prior[1:D], prec = prec_mu[1:D, 1:D])
  }
  
  # Priors for state-specific precision matrices
  for(u in 1:K) {
    Sigma_inv[u, 1:D, 1:D] ~ dwish(R = R_prior[1:D, 1:D], df = alpha_prior)
    Sigma[u, 1:D, 1:D] <- inverse(Sigma_inv[u, 1:D, 1:D])
  }
  
  # Prior for degrees of freedom
  nu_minus_2 ~ dexp(lambda_prior)
  nu <- nu_minus_2 + 2
  
  U[1] ~ dcat(pi[1:K])
  
  for(t in 2:T) {
    U[t] ~ dcat(Pi[U[t-1], 1:K])
  }
  
  # Scaling variables 
  for(t in 1:T) {
    delta[t] ~ dgamma(shape = nu/2, rate = nu/2)
  }
  
  
  for(t in 1:T) {
    for(i in 1:D) {
      for(j in 1:D) {
        Prec_scaled[t, i, j] <- delta[t] * Sigma_inv[U[t], i, j]
      }
    }
    
    Y[t, 1:D] ~ dmnorm(mu[U[t], 1:D], prec = Prec_scaled[t, 1:D, 1:D])
  }
})

# Model setup with hyperparameters

setup_tHMM <- function(Y, K = 4, lambda = 1) {
  
  T <- nrow(Y)
  D <- ncol(Y)
  
  half_alpha <- ceiling((D + 1) / 2) + 1
  alpha <- half_alpha * 2
  
  S <- matrix(0, D, D)
  for(i in 1:D) {
    S[i, i] <- alpha
    if(i < D) {
      for(j in (i+1):D) {
        S[i, j] <- half_alpha
        S[j, i] <- half_alpha
      }
    }
  }
  
  print(S[1:min(5, D), 1:min(5, D)])
  
  R_prior <- solve(S)
  
  # Wilshart matrix (first 5x5)
  
  print(round(R_prior[1:min(5, D), 1:min(5, D)], 4))
  
  # Prior precision for means
  prec_mu <- solve(diag(100, D))
  
  constants <- list(
    K = K,
    D = D,
    T = T,
    alpha = rep(1, K),
    m_prior = rep(0, D),
    prec_mu = prec_mu,
    alpha_prior = alpha,
    R_prior = R_prior,
    lambda_prior = lambda
  )
  
  data <- list(Y = Y)
  
  #K-means initialization
  set.seed(123)
  km <- kmeans(Y, centers = K, nstart = 20)
  U_init <- km$cluster
  
  print(table(U_init))
  
  mu_init <- matrix(0, K, D)
  Sigma_inv_init <- array(0, dim = c(K, D, D))
  
  for(k in 1:K) {
    state_obs <- which(U_init == k)
    if(length(state_obs) > D) {
      mu_init[k, ] <- colMeans(Y[state_obs, , drop = FALSE])
      cov_k <- cov(Y[state_obs, , drop = FALSE])
      
      eigen_vals <- eigen(cov_k)$values
      if(min(eigen_vals) < 0.01) {
        cov_k <- cov_k + diag(0.1, D)
      }
      
      Sigma_inv_init[k, , ] <- solve(cov_k)
    } else {
      mu_init[k, ] <- colMeans(Y)
      Sigma_inv_init[k, , ] <- solve(cov(Y) + diag(0.1, D))
    }
  }
  
  inits <- list(
    pi = rep(1/K, K),
    Pi = matrix(1/K, K, K),
    mu = mu_init,
    Sigma_inv = Sigma_inv_init,
    U = U_init,
    delta = rep(1, T),
    nu_minus_2 = 98  
  )
  
  return(list(
    constants = constants,
    data = data,
    inits = inits,
    hyperparams = list(alpha = alpha, S = S, half_alpha = half_alpha)
  ))
}


# Custom sampler for nu 
sampler_nu_logRW <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    sigma_prop <- 0.5
    calcNodes <- model$getDependencies(target)
    
    # Likelihood of gamma variables
    delta_nodes <- model$expandNodeNames('delta')
  },
  
  run = function() {
    nu_minus_2_current <- model[[target]]
    nu_current <- nu_minus_2_current + 2
    
    theta_current <- log(nu_minus_2_current)
    theta_prop <- rnorm(1, theta_current, sigma_prop)
    nu_minus_2_prop <- exp(theta_prop)
    nu_prop <- nu_minus_2_prop + 2
    
    model[[target]] <<- nu_minus_2_prop
    model_nu <- nu_prop  
    model$calculate('nu')
    
    # Log acceptance ratio
    # Prior ratio
    logPriorRatio <- model$calculateDiff(target)
    
    # Likelihood ratio through delta nodes
    logLikRatio <- model$calculateDiff(delta_nodes)
    
    # Jacobian correction for transformation
    logJacobian <- theta_prop - theta_current
    
    # Total log acceptance ratio
    logMHR <- logPriorRatio + logLikRatio + logJacobian
    
    accept <- decide(logMHR)
    
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, 
           nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, 
           nodes = calcNodes, logProb = TRUE)
      model[[target]] <<- nu_minus_2_current
    }
  },
  
  methods = list(
    reset = function() {}
  )
)



# MCMC with configuration

run_tHMM_MCMC <- function(Y, K = 4, niter = 20000, nburnin = 2000,
                          thin = 1, nchains = 1, sigma_nu = 0.5,
                          id_constraint = 6) {
  
  
  
  model_setup <- setup_tHMM(Y, K)
  
  Rmodel <- nimbleModel(
    code = tHMM_code,
    constants = model_setup$constants,
    data = model_setup$data,
    inits = model_setup$inits
  )
  
  
  Cmodel <- compileNimble(Rmodel)
  
  
  mcmcConf <- configureMCMC(Rmodel, print = TRUE)
  
  mcmcConf$addMonitors(c('Sigma', 'U', 'delta', 'nu', 'pi'))
  
  mcmcConf$removeSamplers('nu_minus_2')
  mcmcConf$addSampler(
    target = 'nu_minus_2',
    type = sampler_nu_logRW,
    control = list(sigma = sigma_nu)
  )
  
  
  Rmcmc <- buildMCMC(mcmcConf)
  
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  
  
  start_time <- Sys.time()
  
  samples <- runMCMC(
    Cmcmc,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    nchains = nchains,
    samplesAsCodaMCMC = TRUE,
    summary = FALSE
  )
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  
  return(list(
    samples = samples,
    model = Cmodel,
    Y = Y,
    K = K,
    hyperparams = model_setup$hyperparams,
    elapsed_time = elapsed
  ))
}


# Label switching

relabel_samples <- function(samples, K, D, constraint_dim = 6) {
  
  
  if(inherits(samples, "mcmc")) {
    samp_mat <- as.matrix(samples)
  } else if(inherits(samples, "mcmc.list")) {
    samp_mat <- as.matrix(samples[[1]])
  } else {
    samp_mat <- samples
  }
  
  n_samples <- nrow(samp_mat)
  
  Sigma_samples <- array(0, dim = c(n_samples, K, D, D))
  
  for(iter in 1:n_samples) {
    for(k in 1:K) {
      for(i in 1:D) {
        for(j in 1:D) {
          col_name <- paste0("Sigma[", k, ", ", i, ", ", j, "]")
          if(col_name %in% colnames(samp_mat)) {
            Sigma_samples[iter, k, i, j] <- samp_mat[iter, col_name]
          }
        }
      }
    }
  }
  
  relabeled_samp <- samp_mat
  n_relabeled <- 0
  
  for(iter in 1:n_samples) {
    variances <- sapply(1:K, function(k) {
      Sigma_samples[iter, k, constraint_dim, constraint_dim]
    })
    
    ord <- order(variances)
    
    if(!all(ord == 1:K)) {
      n_relabeled <- n_relabeled + 1
      
      for(k in 1:K) {
        for(d in 1:D) {
          old_col <- paste0("mu[", ord[k], ", ", d, "]")
          new_col <- paste0("mu[", k, ", ", d, "]")
          if(old_col %in% colnames(samp_mat)) {
            relabeled_samp[iter, new_col] <- samp_mat[iter, old_col]
          }
        }
      }
      
      for(k in 1:K) {
        for(i in 1:D) {
          for(j in 1:D) {
            old_col_sig <- paste0("Sigma[", ord[k], ", ", i, ", ", j, "]")
            new_col_sig <- paste0("Sigma[", k, ", ", i, ", ", j, "]")
            old_col_inv <- paste0("Sigma_inv[", ord[k], ", ", i, ", ", j, "]")
            new_col_inv <- paste0("Sigma_inv[", k, ", ", i, ", ", j, "]")
            
            if(old_col_sig %in% colnames(samp_mat)) {
              relabeled_samp[iter, new_col_sig] <- samp_mat[iter, old_col_sig]
            }
            if(old_col_inv %in% colnames(samp_mat)) {
              relabeled_samp[iter, new_col_inv] <- samp_mat[iter, old_col_inv]
            }
          }
        }
      }
      
      for(i in 1:K) {
        for(j in 1:K) {
          old_col <- paste0("Pi[", ord[i], ", ", ord[j], "]")
          new_col <- paste0("Pi[", i, ", ", j, "]")
          if(old_col %in% colnames(samp_mat)) {
            relabeled_samp[iter, new_col] <- samp_mat[iter, old_col]
          }
        }
      }
      
      for(k in 1:K) {
        old_col <- paste0("pi[", ord[k], "]")
        new_col <- paste0("pi[", k, "]")
        if(old_col %in% colnames(samp_mat)) {
          relabeled_samp[iter, new_col] <- samp_mat[iter, old_col]
        }
      }
      
      for(t in 1:ncol(relabeled_samp)) {
        if(grepl("^U\\[", colnames(relabeled_samp)[t])) {
          old_state <- samp_mat[iter, t]
          new_state <- which(ord == old_state)
          relabeled_samp[iter, t] <- new_state
        }
      }
    }
    
    
  }
  
  cat(sprintf("\n Relabeled %d out of %d iterations (%.1f%%)\n",
              n_relabeled, n_samples, 100*n_relabeled/n_samples))
  
  relabeled_samples <- coda::as.mcmc(relabeled_samp)
  
  return(relabeled_samples)
}



# Posterior estimates extraction
get_posterior_estimates <- function(samples, K, D, var_names = NULL) {
  
  if(inherits(samples, "mcmc") || inherits(samples, "mcmc.list")) {
    samp_mat <- as.matrix(samples)
  } else {
    samp_mat <- samples
  }
  
  estimates <- list()
  
  estimates$mu <- matrix(0, K, D)
  for(k in 1:K) {
    for(d in 1:D) {
      col_name <- paste0("mu[", k, ", ", d, "]")
      if(col_name %in% colnames(samp_mat)) {
        estimates$mu[k, d] <- mean(samp_mat[, col_name])
      }
    }
  }
  
  if(!is.null(var_names)) {
    colnames(estimates$mu) <- var_names
  }
  rownames(estimates$mu) <- paste0("State", 1:K)
  
  estimates$Sigma <- array(0, dim = c(K, D, D))
  
  for(k in 1:K) {
    for(i in 1:D) {
      for(j in 1:D) {
        col_name <- paste0("Sigma[", k, ", ", i, ", ", j, "]")
        if(col_name %in% colnames(samp_mat)) {
          estimates$Sigma[k, i, j] <- mean(samp_mat[, col_name])
        }
      }
    }
    
    if(!is.null(var_names)) {
      dimnames(estimates$Sigma)[[2]] <- var_names
      dimnames(estimates$Sigma)[[3]] <- var_names
    }
  }
  
  estimates$Pi <- matrix(0, K, K)
  for(i in 1:K) {
    for(j in 1:K) {
      col_name <- paste0("Pi[", i, ", ", j, "]")
      if(col_name %in% colnames(samp_mat)) {
        estimates$Pi[i, j] <- mean(samp_mat[, col_name])
      }
    }
  }
  rownames(estimates$Pi) <- paste0("From", 1:K)
  colnames(estimates$Pi) <- paste0("To", 1:K)
  
  if("nu" %in% colnames(samp_mat)) {
    estimates$nu <- mean(samp_mat[, "nu"])
    estimates$nu_sd <- sd(samp_mat[, "nu"])
  }
  
  estimates$pi <- numeric(K)
  for(k in 1:K) {
    col_name <- paste0("pi[", k, "]")
    if(col_name %in% colnames(samp_mat)) {
      estimates$pi[k] <- mean(samp_mat[, col_name])
    }
  }
  names(estimates$pi) <- paste0("State", 1:K)
  
  T <- sum(grepl("^U\\[", colnames(samp_mat)))
  estimates$U <- numeric(T)
  
  for(t in 1:T) {
    u_col <- paste0("U[", t, "]")
    if(u_col %in% colnames(samp_mat)) {
      state_table <- table(samp_mat[, u_col])
      estimates$U[t] <- as.numeric(names(state_table)[which.max(state_table)])
    }
  }
  
  return(estimates)
}



print_results <- function(estimates, K, var_names) {
  
  print(round(estimates$mu, 3))
  
  for(k in 1:K) {
    variances <- diag(estimates$Sigma[k, , ])
    names(variances) <- var_names
    cat(sprintf("State %d:\n", k))
    print(round(variances, 3))
  }
  
  for(k in 1:K) {
    cat(sprintf("\nState %d:\n", k))
    print(round(estimates$Sigma[k, , ], 3))
  }
  
  print(round(estimates$Pi, 3))
  
  cat(sprintf("nu = %.3f (SD = %.3f)\n", estimates$nu, estimates$nu_sd))
  
  print(round(estimates$pi, 3))
  
  print(table(factor(estimates$U, levels = 1:K)))
  
  for(k in 1:K) {
    det_k <- det(estimates$Sigma[k, , ])
    cat(sprintf("State %d: %.3e\n", k, det_k))
  }
  
  for(k in 1:K) {
    cat(sprintf("\nState %d:\n", k))
    D_show <- min(5, length(var_names))
    Sigma_sub <- estimates$Sigma[k, 1:D_show, 1:D_show]
    D_k <- diag(1/sqrt(diag(Sigma_sub)))
    Corr_k <- D_k %*% Sigma_sub %*% D_k
    rownames(Corr_k) <- colnames(Corr_k) <- var_names[1:D_show]
    print(round(Corr_k, 3))
  }
}



forecast_tHMM <- function(mcmc_output, h = 10, n_samples = 1000) {
  
  samples <- mcmc_output$samples
  K <- mcmc_output$K
  Y <- mcmc_output$Y
  D <- ncol(Y)
  T_obs <- nrow(Y)
  
  if(inherits(samples, "mcmc") || inherits(samples, "mcmc.list")) {
    samp_mat <- as.matrix(samples)
  } else {
    samp_mat <- samples
  }
  
  n_mcmc <- nrow(samp_mat)
  if(n_samples > n_mcmc) n_samples <- n_mcmc
  sample_state_obs <- sample(1:n_mcmc, n_samples, replace = FALSE)
  
  Y_forecast <- array(NA, dim = c(n_samples, h, D))
  U_forecast <- matrix(NA, n_samples, h)
  delta_forecast <- matrix(NA, n_samples, h)
  
  
  for(s in 1:n_samples) {
    state_obs <- sample_state_obs[s]
    
    mu <- matrix(0, K, D)
    Pi <- matrix(0, K, K)
    Sigma <- array(0, dim = c(K, D, D))
    
    for(k in 1:K) {
      for(d in 1:D) {
        mu[k, d] <- samp_mat[state_obs, paste0("mu[", k, ", ", d, "]")]
      }
      
      for(i in 1:K) {
        Pi[k, i] <- samp_mat[state_obs, paste0("Pi[", k, ", ", i, "]")]
      }
      
      for(i in 1:D) {
        for(j in 1:D) {
          Sigma[k, i, j] <- samp_mat[state_obs, paste0("Sigma[", k, ", ", i, ", ", j, "]")]
        }
      }
    }
    
    if("nu" %in% colnames(samp_mat)) {
      nu <- samp_mat[state_obs, "nu"]
    } else {
      nu <- samp_mat[state_obs, "nu_minus_2"] + 2
    }
    
    U_current <- samp_mat[state_obs, paste0("U[", T_obs, "]")]
    
    for(t_ahead in 1:h) {
      U_next <- sample(1:K, size = 1, prob = Pi[U_current, ])
      U_forecast[s, t_ahead] <- U_next
      
      delta <- rgamma(1, shape = nu/2, rate = nu/2)
      delta_forecast[s, t_ahead] <- delta
      
      Sigma_scaled <- Sigma[U_next, , ] / delta
      
      eigen_decomp <- eigen(Sigma_scaled, symmetric = TRUE)
      if(min(eigen_decomp$values) < 1e-8) {
        Sigma_scaled <- Sigma_scaled + diag(1e-6, D)
      }
      
      Y_forecast[s, t_ahead, ] <- mvrnorm(n = 1,
                                          mu = mu[U_next, ],
                                          Sigma = Sigma_scaled)
      
      U_current <- U_next
    }
    
  }
  
  
  return(list(
    Y_forecast = Y_forecast,
    U_forecast = U_forecast,
    delta_forecast = delta_forecast,
    h = h,
    n_samples = n_samples,
    D = D,
    variable_names = colnames(Y)
  ))
}

summarize_forecasts <- function(forecast_result) {
  
  h <- forecast_result$h
  D <- forecast_result$D
  Y_forecast <- forecast_result$Y_forecast
  var_names <- forecast_result$variable_names
  
  Y_mean <- apply(Y_forecast, c(2, 3), mean)
  Y_var <- apply(Y_forecast, c(2, 3), var)
  Y_sd <- sqrt(Y_var)
  
  Y_lower <- apply(Y_forecast, c(2, 3), quantile, probs = 0.025)
  Y_upper <- apply(Y_forecast, c(2, 3), quantile, probs = 0.975)
  
  colnames(Y_mean) <- colnames(Y_var) <- colnames(Y_sd) <- var_names
  colnames(Y_lower) <- colnames(Y_upper) <- var_names
  rownames(Y_mean) <- rownames(Y_var) <- rownames(Y_sd) <- paste0("T+", 1:h)
  rownames(Y_lower) <- rownames(Y_upper) <- paste0("T+", 1:h)
  
  U_mode <- apply(forecast_result$U_forecast, 2, function(x) {
    as.numeric(names(which.max(table(x))))
  })
  
  return(list(
    mean = Y_mean,
    variance = Y_var,
    sd = Y_sd,
    lower_95 = Y_lower,
    upper_95 = Y_upper,
    state_mode = U_mode
  ))
}



compute_risk_measures <- function(forecast_result, variable_index = 1, alpha = 0.05) {
  
  h <- forecast_result$h
  Y_forecast <- forecast_result$Y_forecast[, , variable_index]
  var_name <- forecast_result$variable_names[variable_index]
  
  
  VaR <- apply(Y_forecast, 1, quantile, probs = alpha)
  
  ES <- sapply(1:h, function(t) {
    forecasts_t <- Y_forecast[, t]
    mean(forecasts_t[forecasts_t <= VaR[t]])
  })
  
  risk_table <- data.frame(
    Horizon = paste0("T+", 1:h),
    VaR = round(VaR, 3),
    ES = round(ES, 3)
  )
  
  print(risk_table)
  
  return(list(VaR = VaR, ES = ES))
}

plot_state_allocation <- function(Y, estimates, dates = NULL, main_title = NULL) {
  
  T <- nrow(Y)
  D <- ncol(Y)
  
  if (is.null(dates)) {
    Time_vec <- 1:T
  } else {
    if (length(dates) != T) {
      stop("Length of 'dates' must match number of rows in Y.")
    }
    Time_vec <- dates
  }
  
  
  plot_data <- data.frame(
    Time = rep(Time_vec, D),
    Value = as.vector(Y),
    Variable = rep(colnames(Y), each = T),
    State = factor(rep(estimates$U, D))
  )
  
  p <- ggplot(plot_data, aes(x = Time, y = Value, color = State)) +
    geom_point(alpha = 0.5, size = 0.5) +
    facet_wrap(~ Variable, ncol = 1, scales = "free_y") +
    theme_bw() +
    labs(title = if(is.null(main_title)) "Log-Returns with Hidden State Classification" else main_title,
         x = "Time",
         y = "Log-Returns (%)") +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 10))
  
  return(p)
}

plot_forecast_density <- function(forecast_result, variable_index = 1,
                                  horizons = c(1, 5, 10)) {
  
  Y_forecast <- forecast_result$Y_forecast[, , variable_index]
  var_name <- forecast_result$variable_names[variable_index]
  
  plot_data <- data.frame()
  for(h in horizons) {
    if(h <= ncol(Y_forecast)) {
      df_h <- data.frame(
        value = Y_forecast[, h],
        horizon = paste0("T+", h)
      )
      plot_data <- rbind(plot_data, df_h)
    }
  }
  
  p <- ggplot(plot_data, aes(x = value, fill = horizon)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~ horizon, ncol = 1, scales = "free_y") +
    theme_bw() +
    labs(title = paste("Forecast Density for", var_name),
         x = "Log-Return (%)",
         y = "Density") +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Set2")
  
  return(p)
}


K <- 4
niter <- 20000  
nburnin <- 2000  
thin <- 5
sigma_nu <- 0.5
id_constraint <- 6

mcmc_output <- run_tHMM_MCMC(
  Y = Y,
  K = K,
  niter = niter,
  nburnin = nburnin,
  thin = thin,
  nchains = 1,
  sigma_nu = sigma_nu,
  id_constraint = id_constraint
)


relabeled_samples <- relabel_samples(
  samples = mcmc_output$samples,
  K = K,
  D = ncol(Y),
  constraint_dim = id_constraint
)

estimates <- get_posterior_estimates(
  samples = relabeled_samples,
  K = K,
  D = ncol(Y),
  var_names = colnames(Y)
)

print_results(estimates, K, colnames(Y))


ess <- effectiveSize(relabeled_samples)


# mu (first 12)

mu_ess <- ess[grepl("^mu\\[", names(ess))]
mu_ess_matrix <- matrix(mu_ess, nrow = 4, ncol = 7, byrow = TRUE)
rownames(mu_ess_matrix) <- paste0("State ", 1:4)
colnames(mu_ess_matrix) <- c("NVDA", "MSFT", "META", "AAPL", "GOOG", "SP500", "NDX")
print(round(mu_ess_matrix, 2))

#samp_mat <- as.matrix(relabeled_samples)

ess_trimm <- function(samples_mcmc, drop_frac=0.10){
  
  mat <- as.matrix(samples_mcmc) 
  n <- nrow(mat)
  start <- floor(n*drop_frac)+1 
  trimmed <- coda::as.mcmc(mat[start:n, , drop=FALSE])
  coda::effectiveSize(trimmed)
} 

ess_10 <- ess_trimm(relabeled_samples, drop_frac = 0.10) 

ess_25 <- ess_trimm(relabeled_samples, drop_frac = 0.25)

mu_ess_10 <- ess_10[grepl("^mu\\[", names(ess_10))]

mu_ess_25 <- ess_25[grepl("^mu\\[", names(ess_25))]

mu_ess_10_matrix <- matrix(mu_ess_10, nrow = 4, ncol = 7, byrow = TRUE)
mu_ess_25_matrix <- matrix(mu_ess_25, nrow = 4, ncol = 7, byrow = TRUE)

rownames(mu_ess_10_matrix) <- rownames(mu_ess_25_matrix) <- paste0("State ", 1:4)
colnames(mu_ess_10_matrix) <- colnames(mu_ess_25_matrix) <- c("NVDA", "MSFT", "META", "AAPL", "GOOG", "SP500", "NDX")

print(round(mu_ess_10_matrix, 2))

print(round(mu_ess_25_matrix, 2))

traceplot_trimmed_mu <- function(samples_mcmc, drop_frac = 0.10, n_plot = 12, main_title = NULL) {  
  
  mat <- as.matrix(samples_mcmc)
  n <- nrow(mat) 
  start <- floor(n * drop_frac) + 1 
  trimmed <- mat[start:n, , drop = FALSE] 
  
  mu_cols <- grep("^mu\\[", colnames(trimmed), value = TRUE) 
  mu_cols <- head(mu_cols, n_plot) 
  
  op <- par(mfrow = c(4, 3), mar = c(2.5, 2.5, 2, 1), oma = c(0, 0, 2, 0)) 
  
  for (nm in mu_cols) { 
    plot(trimmed[, nm], type = "l", main = nm, xlab = "", ylab = "",
         col = "steelblue") 
  } 
  
  if (is.null(main_title)) {
    main_title <- paste0("Traceplots mu (drop ", drop_frac * 100, "%)") 
  } 
  mtext(main_title, outer = TRUE, cex = 1.1) 
  
  par(op) 
  invisible(trimmed) 
} 

traceplot_trimmed_mu(
  samples_mcmc = relabeled_samples, 
  drop_frac = 0.10, 
  main_title = "Traceplots mu - drop 10%"
) 

traceplot_trimmed_mu(
  samples_mcmc = relabeled_samples, 
  drop_frac = 0.25, 
  main_title = "Traceplots mu - drop 25%"
)


# Pi parameters (first 16)

pi_ess <- ess[grepl("^Pi\\[", names(ess))]
print(head(pi_ess, 16))

cat("\n  nu:", ess["nu"], "\n")

# Sigma parameters (first 12)

sigma_ess <- ess[grepl("^Sigma\\[", names(ess))]
print(head(sigma_ess, 12))

# Plot state allocation 
p_states <- plot_state_allocation(
  Y, 
  estimates,
  dates = dates_clean,
  main_title = "Tech Stock Log-Returns with State Classification"
)
print(p_states)

mcmc_output$samples <- relabeled_samples

forecasts <- forecast_tHMM(mcmc_output, h = 10, n_samples = 500)
forecast_summary <- summarize_forecasts(forecasts)


print(round(forecast_summary$mean, 3))

print(round(forecast_summary$sd, 3))

# Forecast 95% HPR (Lower)
print(round(forecast_summary$lower_95, 3))

# Forecast 95% HPR (Upper)
print(round(forecast_summary$upper_95, 3))

# Risk measures for each stock
for(i in 1:ncol(Y)) {
  risk_measures <- compute_risk_measures(forecasts, variable_index = i, alpha = 0.05)
}

# Plot forecast density for each stock 
p_forecast <- plot_forecast_density(forecasts, variable_index = 1,
                                    horizons = c(1, 5, 10))
print(p_forecast)


forecast_table_stocks <- function(forecast_summary, variable_index, variable_name, h = 10) {
  
  forecast_table <- data.frame(
    Asset = rep(variable_name, h),
    Time = paste0("T+", 1:h),
    E_Y = round(forecast_summary$mean[1:h, variable_index], 3),
    V_Y = round(forecast_summary$variance[1:h, variable_index], 3),
    HPR_lower = round(forecast_summary$lower_95[1:h, variable_index], 3),
    HPR_upper = round(forecast_summary$upper_95[1:h, variable_index], 3)
  )
  
  colnames(forecast_table) <- c("Asset", "Time", 
                                paste0("E[Y_{", variable_index, "T+h}|y_T]"),
                                paste0("V[Y_{", variable_index, "T+h}|y_T]"),
                                "HPR_{0.95} Lower", "HPR_{0.95} Upper")
  
  return(forecast_table)
}

forecasts <- forecast_tHMM(mcmc_output, h = 10, n_samples = 500)
forecast_summary <- summarize_forecasts(forecasts)

stock_names <- colnames(Y)

for(i in 1:ncol(Y)) {
  
  table_i <- forecast_table_stocks(forecast_summary, 
                                   variable_index = i, 
                                   variable_name = stock_names[i],
                                   h = 10)
  print(table_i, row.names = FALSE)
}

save(
  mcmc_output,
  relabeled_samples,
  estimates,
  forecasts,
  forecast_summary,
  Y,
  dates_clean,
  K,
  file = "risultati_20k.RData"
)

load("risultati_20k.RData")


