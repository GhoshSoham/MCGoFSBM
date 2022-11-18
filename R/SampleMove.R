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

  # Output:
  # the graph after one random move
  return(G_sample)
}
