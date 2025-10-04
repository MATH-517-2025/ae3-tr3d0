# nolint start: snake_case_linter, implicit_return_linter
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(shiny)
library(broom)

#set.seed(3009)
# Funziopna M
m <- function(x) {
  sin((x / 3 + 0.1)^-1)
}
# grazie photomath
m_double_prime <- function(x) {
  (2 * (100 * x^2 + 60 * x + 9) * cos(30 / (10 * x + 3))
   - 900 * (x / 3 + 0.1) * sin(30 / (10 * x + 3))) /
    (9 * (10 * x + 3)^2 * (x / 3 + 0.1)^3)
}
#genera campione
generate_sample <- function(n, alpha, beta, sigma2 = 1) {
  X <- rbeta(n, alpha, beta)
  Y <- m(X) + rnorm(n, 0, sqrt(sigma2))
  data.frame(X = X, Y = Y)
}
# un polinomio per ogni blocco, occhio a alfa piccolo e beta grande
fit_block_models <- function(data, N) {
  # Create blocks based on quantiles of X
  breaks <- quantile(data$X, probs = seq(0, 1, length.out = N + 1))
  data$block <- cut(data$X, breaks = breaks, include.lowest = TRUE)

  models <- list()
  block_data <- list()

  for (j in 1:N) {
    block_obs <- data[data$block == levels(data$block)[j], ]
    if (nrow(block_obs) >= 5) {
      # Servono almeno 5 osservazioni per stimare il polinomio 4 p=5
      model <- lm(Y ~ poly(X, 4, raw = TRUE), data = block_obs)
      models[[j]] <- model
      block_data[[j]] <- block_obs
    } else {
      models[[j]] <- NULL
      block_data[[j]] <- NULL
    }
  }
  list(models = models, block_data = block_data, breaks = breaks)
}
# stima theta_22 and sigma^2
estimate_parameters <- function(data, block_results) {
  models <- block_results$models
  block_data <- block_results$block_data
  N <- length(models)

  theta_22_est <- 0
  sigma2_est <- 0
  total_obs <- 0

  for (j in 1:N) {
    if (!is.null(models[[j]])) {
      block_obs <- block_data[[j]]
      n_block <- nrow(block_obs)

      # calcoliamo a mano le derivate, per semplificarci la vita
      # m(x) = b0 + b1*x + b2*x^2 + b3*x^3 + b4*x^4
      # m' (x) = b1 + 2*b2*x + 3*b3*x^2 + 4*b4*x^3
      # m''(x) = 2*b2 + 6*b3*x + 12*b4*x^2
      coefs <- coef(models[[j]])
      if (length(coefs) == 5) {  # Full quartic polynomial
        b2 <- coefs[3]
        b3 <- coefs[4]
        b4 <- coefs[5]

        m_double_prime_est <- 2 * b2 + 6 * b3 * block_obs$X +
          12 * b4 * block_obs$X^2

        theta_22_est <- theta_22_est + sum(m_double_prime_est^2)

        # Calculate residuals for sigma^2
        predictions <- predict(models[[j]], newdata = block_obs)
        residuals <- block_obs$Y - predictions
        sigma2_est <- sigma2_est + sum(residuals^2)

        total_obs <- total_obs + n_block
      }
    }
  }

  # supporto beta [0,1]
  support_length <- 1

  list(
    theta_22 = theta_22_est / total_obs,
    sigma2 = sigma2_est / (total_obs - 5 * N),
    support_length = support_length,
    effective_N = length(compact(models))
  )
}
# optimal bandwidth
calculate_h_amise <- function(n, theta_22, sigma2, support_length) {
  if (is.na(theta_22) || theta_22 <= 0) {
    return(NA)
  }
  n^(-1 / 5) * (35 * sigma2 * support_length / theta_22)^(1 / 5)
}
# Mallow's Cp to find optimal N
find_optimal_N <- function(data, N_max) {
  cp_values <- numeric(N_max)
  rss_values <- numeric(N_max)

  for (N in 1:N_max) {
    block_results <- fit_block_models(data, N)
    param_estimates <- estimate_parameters(data, block_results)

    if (!is.na(param_estimates$sigma2)) {
      rss_values[N] <- param_estimates$sigma2 * (nrow(data) - 5 * N)

      # Cp
      if (N == N_max) {
        rss_Nmax <- rss_values[N_max]
        for (k in 1:N_max) {
          cp_values[k] <- rss_values[k] / (rss_Nmax / (nrow(data) - 5 * N_max))
          - (nrow(data) - 10 * k)
        }
      }
    }
  }

  optimal_N <- which.min(cp_values[1:N_max])
  list(optimal_N = optimal_N, cp_values = cp_values, rss_values = rss_values)
}

# Funzione per studiare l'effetto della dimensione del campione con ripetizioni
study_sample_size_with_reps <- function(n_reps = numero) {
  n_values <- c(100, 200, 500, 1000, 2000, 5000)
  alpha <- 2
  beta <- 2
  sigma2 <- 1
  N_fixed <- 5

  results_list <- list()
  
  for (n in n_values) {
    cat("Processing n =", n, "with", n_reps, "replications\n")
    
    h_amise_values <- numeric(n_reps)
    theta_22_values <- numeric(n_reps)
    sigma2_est_values <- numeric(n_reps)
    
    for (rep in 1:n_reps) {
      set.seed(3009 + n * 1000 + rep)  # Seed variabile per ogni combinazione
      data <- generate_sample(n, alpha, beta, sigma2)
      block_results <- fit_block_models(data, N_fixed)
      param_estimates <- estimate_parameters(data, block_results)
      
      h_amise <- calculate_h_amise(n, param_estimates$theta_22,
                                   param_estimates$sigma2,
                                   param_estimates$support_length)
      
      h_amise_values[rep] <- h_amise
      theta_22_values[rep] <- param_estimates$theta_22
      sigma2_est_values[rep] <- param_estimates$sigma2
    }
    
    # Calcola medie e deviazioni standard
    results_list[[as.character(n)]] <- data.frame(
      n = n,
      h_amise_mean = mean(h_amise_values, na.rm = TRUE),
      h_amise_sd = sd(h_amise_values, na.rm = TRUE),
      theta_22_mean = mean(theta_22_values, na.rm = TRUE),
      theta_22_sd = sd(theta_22_values, na.rm = TRUE),
      sigma2_est_mean = mean(sigma2_est_values, na.rm = TRUE),
      sigma2_est_sd = sd(sigma2_est_values, na.rm = TRUE),
      n_reps = sum(!is.na(h_amise_values))
    )
  }
  
  results <- do.call(rbind, results_list)
  
  # Plot risultati
  p1 <- ggplot(results, aes(x = n, y = h_amise_mean)) +
    geom_point(size = 3) +
    geom_line() +
    geom_errorbar(aes(ymin = h_amise_mean - h_amise_sd, 
                      ymax = h_amise_mean + h_amise_sd), 
                  width = 0.1) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = paste("Effect of Sample Size on Optimal Bandwidth (", n_reps, "replications)"),
         x = "Sample Size (n)", y = "h_AMISE") +
    theme_minimal()
  
  print(p1)
  ggsave("sample_size_with_reps.png", p1, width = 12, height = 8, dpi = 300)
  
  return(results)
}

# Funzione per studiare l'effetto della dimensione dei blocchi con ripetizioni
study_block_size_with_reps <- function(n_reps = numero) {
  n <- 1000
  alpha <- 2
  beta <- 2
  sigma2 <- 1
  
  # Trova N massimo
  N_max <- max(min(floor(n/20), 5), 1)
  N_values <- 1:N_max
  
  results_list <- list()
  
  for (N in N_values) {
    cat("Processing N =", N, "with", n_reps, "replications\n")
    
    h_amise_values <- numeric(n_reps)
    theta_22_values <- numeric(n_reps)
    sigma2_est_values <- numeric(n_reps)
    optimal_N_values <- numeric(n_reps)
    cp_values <- numeric(n_reps)
    
    for (rep in 1:n_reps) {
      set.seed(3009 + N * 1000 + rep)  # Seed variabile per ogni combinazione
      data <- generate_sample(n, alpha, beta, sigma2)
      block_results <- fit_block_models(data, N)
      param_estimates <- estimate_parameters(data, block_results)
      
      h_amise <- calculate_h_amise(n, param_estimates$theta_22, 
                                  param_estimates$sigma2, param_estimates$support_length)
      
      h_amise_values[rep] <- h_amise
      theta_22_values[rep] <- param_estimates$theta_22
      sigma2_est_values[rep] <- param_estimates$sigma2
      
      # Trova N ottimale per questo campione
      cp_result <- find_optimal_N(data, N_max)
      optimal_N_values[rep] <- cp_result$optimal_N
      cp_values[rep] <- cp_result$cp_values[N]
    }
    
    results_list[[as.character(N)]] <- data.frame(
      N = N,
      h_amise_mean = mean(h_amise_values, na.rm = TRUE),
      h_amise_sd = sd(h_amise_values, na.rm = TRUE),
      theta_22_mean = mean(theta_22_values, na.rm = TRUE),
      theta_22_sd = sd(theta_22_values, na.rm = TRUE),
      sigma2_est_mean = mean(sigma2_est_values, na.rm = TRUE),
      sigma2_est_sd = sd(sigma2_est_values, na.rm = TRUE),
      optimal_N_mean = mean(optimal_N_values, na.rm = TRUE),
      Cp_mean = mean(cp_values, na.rm = TRUE),
      Cp_sd = sd(cp_values, na.rm = TRUE),
      n_reps = sum(!is.na(h_amise_values))
    )
  }
  
  results <- do.call(rbind, results_list)
  
  # Trova N ottimale globale (media tra tutte le repliche)
  global_optimal_N <- round(mean(results$optimal_N_mean))
  
  # Plot risultati
  p1 <- ggplot(results, aes(x = N, y = h_amise_mean)) +
    geom_point(size = 3) +
    geom_line() +
    geom_errorbar(aes(ymin = h_amise_mean - h_amise_sd, 
                      ymax = h_amise_mean + h_amise_sd), # non è un CI, ma soltanto una dispersione di +-1 sd attorno ala media
                  width = 0.1) +
    geom_vline(xintercept = global_optimal_N, linetype = "dashed", color = "red") +
    annotate("text", x = global_optimal_N, y = max(results$h_amise_mean), 
             label = paste("Optimal N =", global_optimal_N), vjust = -1, color = "red") +
    labs(title = paste("Effect of Block Size on Optimal Bandwidth (", n_reps, "replications)"),
         x = "Number of Blocks (N)", y = "h_AMISE") +
    theme_minimal()
  
  print(p1)
  ggsave("block_size_with_reps.png", p1, width = 12, height = 8, dpi = 300)
    cp_data <- data.frame(N = 1:N_max, Cp = cp_result$cp_values[1:N_max])
  p2 <- ggplot(results, aes(x = N, y = Cp_mean)) +
    geom_point(size = 3) +
    geom_line() +
    geom_errorbar(aes(ymin = Cp_mean - Cp_sd, ymax = Cp_mean + Cp_sd), width = 0.1) +
    geom_vline(xintercept = global_optimal_N, linetype = "dashed", color = "red") +
    labs(title = paste("Mallow's Cp for Block Size Selection (", n_reps, "replications)"),
         x = "Number of Blocks (N)", y = "Cp") +
    theme_minimal()
  
  print(p2)
  ggsave("block_size_cp_with_reps.png", p2, width = 12, height = 8, dpi = 300)

  return(list(results = results, global_optimal_N = global_optimal_N))
}

# Funzione per studiare l'effetto dei parametri Beta con ripetizioni
study_beta_params_with_reps <- function(n_reps = numero) {
  n <- 1000
  sigma2 <- 1
  N_fixed <- 5
  
  param_combinations <- expand.grid(
    alpha = c(0.5, 1, 2, 3, 4, 5), 
    beta = c(0.5, 1, 2, 3, 4, 5)
  )
  
  results_list <- list()
  
  for (i in 1:nrow(param_combinations)) {
    alpha <- param_combinations$alpha[i]
    beta <- param_combinations$beta[i]
    
    cat("Processing alpha =", alpha, "beta =", beta, "with", n_reps, "replications\n")
    
    h_amise_values <- numeric(n_reps)
    theta_22_values <- numeric(n_reps)
    sigma2_est_values <- numeric(n_reps)
    
    for (rep in 1:n_reps) {
      set.seed(3009 + alpha * 100 + beta * 10 + rep)  # Seed variabile per ogni combinazione
      data <- generate_sample(n, alpha, beta, sigma2)
      block_results <- fit_block_models(data, N_fixed)
      param_estimates <- estimate_parameters(data, block_results)
      
      h_amise <- calculate_h_amise(n, param_estimates$theta_22, 
                                  param_estimates$sigma2, param_estimates$support_length)
      
      h_amise_values[rep] <- h_amise
      theta_22_values[rep] <- param_estimates$theta_22
      sigma2_est_values[rep] <- param_estimates$sigma2
    }
    
    results_list[[i]] <- data.frame(
      alpha = alpha,
      beta = beta,
      h_amise_mean = mean(h_amise_values, na.rm = TRUE),
      h_amise_sd = sd(h_amise_values, na.rm = TRUE),
      theta_22_mean = mean(theta_22_values, na.rm = TRUE),
      theta_22_sd = sd(theta_22_values, na.rm = TRUE),
      sigma2_est_mean = mean(sigma2_est_values, na.rm = TRUE),
      sigma2_est_sd = sd(sigma2_est_values, na.rm = TRUE),
      n_reps = sum(!is.na(h_amise_values))
    )
  }
  
  results <- do.call(rbind, results_list)
  
  # Plot risultati
  p1 <- ggplot(results, aes(x = factor(alpha), y = factor(beta), fill = h_amise_mean)) +
    geom_tile() +
    geom_text(aes(label = round(h_amise_mean, 3)), color = "white", size = 3) +
    scale_fill_viridis_c() +
    labs(title = paste("Effect of Beta Distribution Parameters on Optimal Bandwidth (", n_reps, "replications)"),
         x = "Alpha", y = "Beta", fill = "h_AMISE") +
    theme_minimal()
  
  print(p1)
  ggsave("beta_params_with_reps.png", p1, width = 12, height = 8, dpi = 300)
  
  return(results)
}

# Esegui le analisi con ripetizioni
sample_size_results_reps <- study_sample_size_with_reps(n_reps = 200)
block_size_results_reps <- study_block_size_with_reps(n_reps = 200)
beta_params_results_reps <- study_beta_params_with_reps(n_reps = 200)

# nolint end