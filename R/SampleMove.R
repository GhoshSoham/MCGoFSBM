# Sample a graph in the same fiber
# The true model is ERSBM

sample_a_move <- function(C, G_current) {
  # Input:
  # G_current: igraph object and undirected graph
  # C: vector of block assignment

  # Getting block information
  num_blocks <- length(unique(C))

  # Getting edge information
  num_edges <- length(E(G_current))
  all_edges <- get.edgelist(G_current)

  # Getting edge information of the complement graph
  G_comp <- graph.complementer(G_current, loops = FALSE)
  comp_edges <- get.edgelist(G_comp)

  # Determine whether the move is inter-block or intra-block
  type <- sample(2, size = 1)

  if (type == 1) {
    # Adding and deleting edges from same block
    # Sample a block
    s <- sample(num_blocks, size = 1)

    # Find vertexes within the block for both the graph and complement graph
    to_delete <- all_edges[((C[as.numeric(all_edges[, 1])] == s) * (C[as.numeric(all_edges[, 2])] == s)) > 0, ]
    to_add <- comp_edges[((C[as.numeric(comp_edges[, 1])] == s) * (C[as.numeric(comp_edges[, 2])] == s)) > 0, ]

    # Check whether the sampled block and complement has at least one edge
    if ((length(to_delete) > 0) * (length(to_add) > 0)) {
      # Sample an edge to add from complement graph and delete from the graph
      delete_edge <- sample(nrow(to_delete), 1)
      add_edge <- sample(nrow(to_add), 1)

      # Add and delete an edge in the same block
      G_sample <- G_current %>%
        add_edges(to_add[add_edge, ]) %>%
        delete_edges(paste0(to_delete[delete_edge, 1], "|", to_delete[delete_edge, 2]))

    } else {
      G_sample <- G_current
    }
  } else if (type == 2) {
    # Adding and deleting edges between two different blocks

  }

  # Output:
  # the graph after one random move
  return(G_sample)
}
