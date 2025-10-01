# nolint start: snake_case_linter, implicit_return_linter
source("bandwidth_selection.R")
# Shiny App
bandwidth_app <- function() {
  ui <- fluidPage(
    titlePanel("Global Bandwidth Selection Simulation"),
    
    sidebarLayout(
      sidebarPanel(
        sliderInput("n", "Sample Size (n):", min = 100, max = 5000, value = 1000, step = 100),
        sliderInput("alpha", "Beta Alpha:", min = 0.5, max = 10, value = 2, step = 0.25),
        sliderInput("beta", "Beta Beta:", min = 0.5, max = 10, value = 2, step = 0.25),
        sliderInput("sigma2", "Error Variance:", min = 0.1, max = 5, value = 1, step = 0.1),
        sliderInput("N", "Number of Blocks (N):", min = 1, max = 10, value = 5, step = 1),
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
      
      # Generate theoretical curve for comparison
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
      N_max <- max(min(floor(input$n/20), 10), 2)
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