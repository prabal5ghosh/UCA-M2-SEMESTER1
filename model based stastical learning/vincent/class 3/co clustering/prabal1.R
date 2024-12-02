simulate_data <- function(param, n_rows, n_cols) {
  # Determine row and column groups
  row_groups <- sample(1:length(param$prop.r), n_rows, replace = TRUE, prob = param$prop.r)
  col_groups <- sample(1:length(param$prop.c), n_cols, replace = TRUE, prob = param$prop.c)
  
  # Initialize the data matrix
  data <- matrix(0, n_rows, n_cols)
  
  # Fill the matrix based on block parameters
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      row_group <- row_groups[i]
      col_group <- col_groups[j]
      block_param <- param$a[[row_group]][[col_group]]
      data[i, j] <- rnorm(1, mean = block_param$mu, sd = sqrt(block_param$var))
    }
  }
  
  return(list(data = data, row_groups = row_groups, col_groups = col_groups))
}

# Simulation with given parameters
set.seed(123) # For reproducibility
sim_result <- simulate_data(param, n_rows = 100, n_cols = 100)
sim_data <- sim_result$data

# Visualize the simulated data
heatmap(sim_data, Rowv = NA, Colv = NA, scale = "none", col = terrain.colors(256))

