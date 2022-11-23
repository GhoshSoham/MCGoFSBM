# The true model is ERSBM
# Chi-Square goodness of fit test statistic computation using the observed graph and MLE

graphchi <- function(G, C, p_mle) {
  # Input: G: igraph object which is an undirected graph and has no self loop
  #        C: numeric vector of size n of block assignment; from 1 to k
  #        p_mle: k*k matrix of MLE table

  # Getting graph information
  n <- length(C)
  k <- length(unique(C))
  A <- igraph::get.adjacency(G, type = "both")

  # Observed no of neighbors of a element in each block
  degseq <- matrix(0, nrow = n, ncol = k)

  # Expected no of neighbors of a element in each block
  exp_mat <- matrix(0, nrow = n, ncol = k)
  for (inode in 1:n) {
    row <- as.numeric(A[inode, ])
    for (iblock in 1:k) {
      degseq[inode, iblock] <- sum(row[C == iblock])
      ni <- sum(C == iblock) # no of nodes in ith block
      exp_mat[inode, iblock] <- ni * p_mle[iblock, C[inode]]
    }
  }

  # Output:
  # the value for chi-sq statistics
  return(sum((degseq[exp_mat != 0] - exp_mat[exp_mat != 0])^2 / exp_mat[exp_mat != 0]))
}
