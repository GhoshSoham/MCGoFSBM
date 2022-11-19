# Sample a graph in the same fiber
# The true model is ERSBM

sample_a_move <- function(C, G_current) {
  # Input:
  # G_current: igraph object which is an undirected graph and has no self loop
  # C: numeric vector of size n of block assignment; from 1 to k

  # Getting block information
  num_blocks <- length(unique(C))

  # Getting edge information
  num_edges <- length(igraph::E(G_current))
  all_edges <- igraph::get.edgelist(G_current)

  # Getting edge information of the complement graph
  G_comp <- igraph::graph.complementer(G_current, loops = FALSE)
  comp_edges <- igraph::get.edgelist(G_comp)

  # Determine whether the move is inter-block or intra-block
  type <- sample.int(2, size = 1)

  if (type == 1) {
    # Adding and deleting edges from same block
    # Sample a block
    s <- sample.int(num_blocks, size = 1)

    # Find edges within the sampled block for both the graph and complement graph
    to_delete <- all_edges[((C[as.numeric(all_edges[, 1])] == s) * (C[as.numeric(all_edges[, 2])] == s)) > 0, ]
    to_add <- comp_edges[((C[as.numeric(comp_edges[, 1])] == s) * (C[as.numeric(comp_edges[, 2])] == s)) > 0, ]

    # Check whether the sampled block has at least one edge
    if ((length(to_delete) > 0) * (length(to_add) > 0)) {
      # Sample an edge to add from complement graph and delete from the graph
      delete_edge <- sample.int(nrow(to_delete), 1)
      add_edge <- sample.int(nrow(to_add), 1)

      # Add and delete an edge in the same block
      G_sample <- igraph::add_edges(G_current, to_add[add_edge, ])
      G_sample <- igraph::delete_edges(G_sample, paste0(to_delete[delete_edge, 1], "|", to_delete[delete_edge, 2]))
    } else {
      G_sample <- G_current
    }
  } else if (type == 2) {
    # Adding and deleting edges between two different blocks
    # Sample two different blocks
    two_blocks <- sample.int(num_blocks, size = 2, replace = FALSE)
    s <- two_blocks[1]
    t <- two_blocks[2]

    # Find edges between two fixed blocks for both the graph and complement graph
    inter <- all_edges[((C[as.numeric(all_edges[, 1])] == s) * (C[as.numeric(all_edges[, 2])] == t)) + ((C[as.numeric(all_edges[, 1])] == t) * (C[as.numeric(all_edges[, 2])] == s)) > 0, ]
    comp_inter <- comp_edges[((C[as.numeric(comp_edges[, 1])] == s) * (C[as.numeric(comp_edges[, 2])] == t)) + ((C[as.numeric(comp_edges[, 1])] == t) * (C[as.numeric(comp_edges[, 2])] == s)) > 0, ]

    # Check whether the sampled blocks have at least one edge
    if ((length(inter) > 0) * (length(comp_inter) > 0)) {
      # Sample an edge to add from complement graph and delete from the graph
      delete_edge <- sample.int(nrow(inter), 1)
      add_edge <- sample.int(nrow(comp_inter), 1)

      # Add and delete an edge between two fixed block
      G_sample <- igraph::add_edges(G_current, comp_inter[add_edge, ])
      G_sample <- igraph::delete_edges(G_sample, paste0(inter[delete_edge, 1], "|", inter[delete_edge, 2]))
    } else {
      G_sample <- G_current
    }
  }

  # Output:
  # the graph after one random move
  return(G_sample)
}
