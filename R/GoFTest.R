# The true model is ERSBM
# Chi-Square goodness of fit test statistic computation using the observed graph and MLE

goftest <- function(G_obs, C, len_chain) {
  # Input: G: igraph object of a graph, undirected graph and no self loop
  #        C: numeric vector of size n of block assignment; from 1 to k
  #        len_chain: number of times the chain will move and generate new graph

  if(class(G_obs) != "igraph"){
    stop("The input G_obs should be an igraph object")
  }

  if(igraph::is.weighted(G_obs)){
    stop("The input graph G_obs shouldn't be a weighted one")
  }

  if(igraph::is_directed(G_obs)){
    stop("The input graph G_obs shouldn't be a directed one")
  }

  if(igraph::any_loop(G_obs)){
    stop("The graph G_obs shouldn't have any self loop")
  }


  # Output:
  # p-value
  return(chi_seq)
}
