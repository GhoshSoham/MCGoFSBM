# The true model is ERSBM
# Estimate the corresponding parameter values i.e,
# q_i,j (the probability of edges between block i and j)


get_mle <- function(G, C) {
  # Input:
  # G_current: G_obs igraph object which is an undirected graph and has no self loop
  # C: numeric vector of size n of block assignment; from 1 to k

  # Getting graph information
  k <- length(unique(C)) # no of block
  n <- length(C) # no of nodes
  A <- igraph::get.adjacency(G, type = "both")

  # Calculate total observed edges between block and within block
  table_obs <- matrix(0, nrow = k, ncol = k)
  for (i in 1:k) {
    for (j in 1:k) {
      table_obs[i, j] <- sum(A[C == i, C == j])
    }
  }

  # Calculate total possible edges between and within block based on nodes information
  tot_edge <- tcrossprod(as.vector(table(C))) - diag(as.vector(table(C)))

  # Estimated value for q_i,j in matrix form
  mle_mat <- table_obs / tot_edge

  # Output:
  # the estimated q_i,j matrix (k by k)
  return(mle_mat)
}
