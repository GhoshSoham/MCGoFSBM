# Sample a graph in the same fiber
# The true model is ERSBM

sample_a_move_cpp <- function(C, A) {
  # Input:
  # A: a n by n binary symmetric adjacency matrix representing a undirected graph
  # C: numeric vector of size n of block assignment; from 1 to k

  # Getting block information
  num_blocks <- length(unique(C))
  n <- length(C)

  # Getting edge information of the graph
  num_edges <- sum(A) / 2
  all_edges <- get_edgelist(A, num_edges, n)

  # Getting edge information of the complement graph
  num_edges_comp <- (n * (n - 1)) / 2 - num_edges
  comp_edges <- get_edgelist_comp(A, num_edges_comp, n)

  # Determine whether the move is inter-block or intra-block
  type <- sample.int(2, size = 1)


  # Output:
  # the graph after one random move
  return(A)
}
