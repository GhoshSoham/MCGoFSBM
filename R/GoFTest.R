# The true model is ERSBM
# Testing whether the observed graph follows ERSBM model using chi square goodness of fit test
# for network data

#'
#' @title Monte Carlo Goodness of Fit tests for Stochastic Block Models
#'
#' @description `goftest` performs chi square goodness of fit test for network data considering the model as ERSBM
#'
#' @param A a n by n binary symmetric adjacency matrix representing a undirected graph where n is the no nodes in the graph
#' @param K a numeric scalar representing no of blocks
#' @param numGraphs number of graphs will be sampled; default value is 100
#'
#' @return A list with the elements
#' \item{statistic}{Value of the chi-square test statistic}
#' \item{p.value}{the p-value for the test.}
#'
#' @export
#'
#' @examples
goftest <- function(A, K, numGraphs = 100) {
  # Some compatibility checks and error message
  # Check whether the input A is a matrix
  if (!is.matrix(A)) {
    stop("A should be an adjacency matrix of the graph.")
  }

  # Check whether the graph corresponding to A is an undirected
  if (!isSymmetric.matrix(A)) {
    stop("A should be a square symmetric matrix.")
  }

  # Check whether the graph corresponding to A is an unweighted
  if (!all(A %in% c(0, 1))) {
    stop("A can only contain 0's and 1's.")
  }

  # Check whether the graph corresponding to A has no self loops
  if (!all(diag(A) == 0)) {
    stop("All the diagonal entries of A should be 0")
  }

  # Getting block assignment for each node from adjacency matrix
  community <- randnet::reg.SP(A, K)
  C <- community$cluster

  # Getting an igraph (an undirected and unweighted) object from the input adjacency matrix
  G_obs <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = NULL)

  # Calculate estimate of the parameter from observed graph
  # which will remain same after generating a new graph on same fiber
  p_mle <- get_mle(G_obs, C)

  # It will store GoF test statistic on graphs
  chi_seq <- rep(0, numGraphs)
  G <- G_obs

  # Storing the first entry of chi_seq as test-stat on observed graph
  chi_seq[1] <- graphchi(G, C, p_mle)
  for (i in 2:numGraphs) {
    # Sampling a new graph
    G_current <- sample_a_move(C, G)

    # Computing GoF test statistic on new sampled graph
    chi_seq[i] <- graphchi(G_current, C, p_mle)
    G <- G_current
  }

  # pvalue i.e, proportion of sampled grpahs has larger GoF statistic than observed one
  pvalue <- mean(chi_seq > chi_seq[1])

  # Output:
  # chi_seq: sequence of chi square test statistics on the sampled graphs
  # pvalue: estimated p-value when true model is ERSBM
  return(list(statistic = chi_seq[1], p.value = pvalue))
}
