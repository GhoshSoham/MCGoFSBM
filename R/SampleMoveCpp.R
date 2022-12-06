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

  if (type == 1) {
    # Adding and deleting edges from same block
    # Sample a block
    s <- sample.int(num_blocks, size = 1)

    # Find edges within the sampled block for both the graph and complement graph
    to_delete <- all_edges[((C[as.numeric(all_edges[, 1])] == s) * (C[as.numeric(all_edges[, 2])] == s)) > 0, ]
    to_add <- comp_edges[((C[as.numeric(comp_edges[, 1])] == s) * (C[as.numeric(comp_edges[, 2])] == s)) > 0, ]

    # Getting dimension information
    to_delete_l <- length(to_delete)
    to_add_l <- length(to_add)

    # Check whether the sampled block has at least one edge
    if ((to_delete_l > 0) * (to_add_l > 0)) {
      # Sample an edge to add from complement graph and delete from the graph
      delete_edge <- sample.int(to_delete_l / 2, 1)
      add_edge <- sample.int(to_add_l / 2, 1)

      # If the sampled block has only one edge in the complement graph, then that will be added,
      # otherwise it will add one edge by sampling randomly from complement graph
      if (length(to_add) == 2) {
        A[to_add[1], to_add[2]] <- 1
        A[to_add[2], to_add[1]] <- 1
      } else {
        A[to_add[add_edge, ][1], to_add[add_edge, ][2]] <- 1
        A[to_add[add_edge, ][2], to_add[add_edge, ][1]] <- 1
      }

      # If the sampled block has only one edge in the graph, then that will be deleted
      # otherwise it will delete one edge by sampling randomly from graph
      if (to_delete_l == 2) {
        A[to_delete[1], to_delete[2]] <- 0
        A[to_delete[2], to_delete[1]] <- 0
      } else {
        A[to_delete[delete_edge, ][1], to_delete[delete_edge, ][2]] <- 0
        A[to_delete[delete_edge, ][2], to_delete[delete_edge, ][1]] <- 0
      }
    }
  } else if (type == 2) {

  }

  # Output:
  # the graph after one random move
  return(A)
}
