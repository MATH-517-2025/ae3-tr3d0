# nolint start: snake_case_linter, implicit_return_linter
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(shiny)
library(broom)

# Funzione M
m <- function(x) {
  sin((x / 3 + 0.1)^-1)
}

# Genera campione
generate_sample <- function(n, alpha, beta, sigma2 = 1) {
  X <- rbeta(n, alpha, beta)
  Y <- m(X) + rnorm(n, 0, sqrt(sigma2))
  data.frame(X = X, Y = Y)
}

# Fit polinomio per ogni blocco
fit_block_models <- function(data, N) {
  breaks <- quantile(data$X, probs = seq(0, 1, length.out = N + 1))
  data$block <- cut(data$X, breaks = breaks, include.lowest = TRUE)

  models <- list()
  block_data <- list()

  for (j in 1:N) {
    block_obs <- data[data$block == levels(data$block)[j], ]
    if (nrow(block_obs) >= 5) {
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

# Stima theta_22 e sigma^2
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

      coefs <- coef(models[[j]])
      if (length(coefs) == 5) {
        b2 <- coefs[3]
        b3 <- coefs[4]
        b4 <- coefs[5]

        m_double_prime_est <- 2 * b2 + 6 * b3 * block_obs$X + 12 * b4 * block_obs$X^2
        theta_22_est <- theta_22_est + sum(m_double_prime_est^2)

        predictions <- predict(models[[j]], newdata = block_obs)
        residuals <- block_obs$Y - predictions
        sigma2_est <- sigma2_est + sum(residuals^2)

        total_obs <- total_obs + n_block
      }
    }
  }

  list(
    theta_22 = theta_22_est / total_obs,
    sigma2 = sigma2_est / (total_obs - 5 * N),
    support_length = 1,
    effective_N = length(compact(models))
  )
}

# Optimal bandwidth
calculate_h_amise <- function(n, theta_22, sigma2, support_length) {
  if (is.na(theta_22) || theta_22 <= 0) {
    return(NA)
  }
  n^(-1 / 5) * (35 * sigma2 * support_length / theta_22)^(1 / 5)
}

# Mallow's Cp per trovare N ottimale
find_optimal_N <- function(data, N_max) {
  cp_values <- numeric(N_max)
  rss_values <- numeric(N_max)

  for (N in 1:N_max) {
    block_results <- fit_block_models(data, N)
    param_estimates <- estimate_parameters(data, block_results)

    if (!is.na(param_estimates$sigma2)) {
      rss_values[N] <- param_estimates$sigma2 * (nrow(data) - 5 * N)

      if (N == N_max) {
        rss_Nmax <- rss_values[N_max]
        for (k in 1:N_max) {
          cp_values[k] <- rss_values[k] / (rss_Nmax / (nrow(data) - 5 * N_max)) - (nrow(data) - 10 * k)
        }
      }
    }
  }

  optimal_N <- which.min(cp_values[1:N_max])
  list(optimal_N = optimal_N, cp_values = cp_values, rss_values = rss_values)
}

# Shiny App
bandwidth_app <- function() {
  ui <- fluidPage(
    titlePanel("Global Bandwidth Selection Simulation"),
    
    sidebarLayout(
      sidebarPanel(
        sliderInput("n", "Sample Size (n):", min = 100, max = 5000, value = 1000, step = 100),
        sliderInput("alpha", "Beta Alpha:", min = 0.5, max = 5, value = 2, step = 0.25),
        sliderInput("beta", "Beta Beta:", min = 0.5, max = 5, value = 2, step = 0.25),
        sliderInput("sigma2", "Error Variance:", min = 0.1, max = 5, value = 1, step = 0.1),
        sliderInput("N", "Number of Blocks (N):", min = 1, max = 5, value = 5, step = 1),
        actionButton("run", "Run Simulation")
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel("Data Overview", 
                   plotOutput("data_plot"),
                   plotOutput("density_plot")),
          tabPanel("Bandwidth Results", 
                   verbatimTextOutput("results_text"),
                   plotOutput("bandwidth_plot")),
          tabPanel("Block Analysis",
                   plotOutput("block_plot"),
                   plotOutput("cp_plot"))
        )
      )
    )
  )
  
  server <- function(input, output) {
    
    simulation_data <- eventReactive(input$run, {
      generate_sample(input$n, input$alpha, input$beta, input$sigma2)
    })
    
    block_results <- eventReactive(input$run, {
      fit_block_models(simulation_data(), input$N)
    })
    
    param_estimates <- eventReactive(input$run, {
      estimate_parameters(simulation_data(), block_results())
    })
    
    output$data_plot <- renderPlot({
      data <- simulation_data()
      ggplot(data, aes(x = X, y = Y)) +
        geom_point(alpha = 0.5) +
        stat_function(fun = m, color = "red", lwd = 1) +
        labs(title = "Simulated Data with True Regression Function",
             x = "X", y = "Y") +
        theme_minimal()
    })
    
    output$density_plot <- renderPlot({
      data <- simulation_data()
      ggplot(data, aes(x = X)) +
        geom_density(fill = "blue", alpha = 0.5) +
        labs(title = "Distribution of Covariate X",
             x = "X", y = "Density") +
        theme_minimal()
    })
    
    output$results_text <- renderPrint({
      estimates <- param_estimates()
      h_amise <- calculate_h_amise(input$n, estimates$theta_22, 
                                  estimates$sigma2, estimates$support_length)
      
      cat("Simulation Results:\n")
      cat("===================\n")
      cat("Sample size (n):", input$n, "\n")
      cat("Beta parameters (α, β):", input$alpha, ",", input$beta, "\n")
      cat("Number of blocks (N):", input$N, "\n")
      cat("Effective blocks:", estimates$effective_N, "\n")
      cat("\nParameter Estimates:\n")
      cat("θ₂₂:", round(estimates$theta_22, 4), "\n")
      cat("σ²:", round(estimates$sigma2, 4), "\n")
      cat("Support length:", estimates$support_length, "\n")
      cat("\nOptimal Bandwidth:\n")
      cat("h_AMISE:", round(h_amise, 4), "\n")
    })
    
    output$bandwidth_plot <- renderPlot({
      estimates <- param_estimates()
      h_amise <- calculate_h_amise(input$n, estimates$theta_22, 
                                  estimates$sigma2, estimates$support_length)
      
      n_theoretical <- seq(100, 5000, length.out = 100)
      h_theoretical <- n_theoretical^(-1/5) * (35 * estimates$sigma2 * 1 / estimates$theta_22)^(1/5)
      
      theoretical_data <- data.frame(n = n_theoretical, h = h_theoretical)
      
      ggplot(theoretical_data, aes(x = n, y = h)) +
        geom_line(color = "blue") +
        geom_point(data = data.frame(n = input$n, h = h_amise), 
                  aes(x = n, y = h), color = "red", lwd = 3) +
        scale_x_log10() +
        scale_y_log10() +
        labs(title = "Theoretical h_AMISE vs Sample Size",
             x = "Sample Size (n)", y = "h_AMISE",
             subtitle = paste("Current estimate: h =", round(h_amise, 4))) +
        theme_minimal()
    })
    
    output$block_plot <- renderPlot({
      data <- simulation_data()
      block_res <- block_results()
      
      plot_data <- data
      plot_data$block <- cut(plot_data$X, breaks = block_res$breaks, include.lowest = TRUE)
      
      ggplot(plot_data, aes(x = X, y = Y, color = block)) +
        geom_point(alpha = 0.6) +
        labs(title = "Data Partitioned into Blocks",
             x = "X", y = "Y", color = "Block") +
        theme_minimal() +
        theme(legend.position = "bottom")
    })
    
    output$cp_plot <- renderPlot({
      data <- simulation_data()
      N_max <- max(min(floor(input$n/20), 5), 2)
      cp_result <- find_optimal_N(data, N_max)
      
      cp_data <- data.frame(N = 1:N_max, Cp = cp_result$cp_values[1:N_max])
      
      ggplot(cp_data, aes(x = N, y = Cp)) +
        geom_point(size = 3) +
        geom_line() +
        geom_vline(xintercept = cp_result$optimal_N, linetype = "dashed", color = "red") +
        labs(title = "Mallow's Cp for Optimal Block Size Selection",
             x = "Number of Blocks (N)", y = "Cp",
             subtitle = paste("Optimal N =", cp_result$optimal_N)) +
        theme_minimal()
    })
  }
  
  shinyApp(ui, server)
}

# Uncomment to run the Shiny app
bandwidth_app()

# nolint end