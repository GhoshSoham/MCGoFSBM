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
#' @param C numeric vector of size n of block assignment of each node; from 1 to K
#' @param numGraphs number of graphs will be sampled; default value is 100
#'
#' @return A list with the elements
#' \item{statistic}{Value of the chi-square test statistic}
#' \item{p.value}{the p-value for the test.}
#'
#' @export
#'
#' @examples
goftest <- function(A, K, C, numGraphs = 100) {
  # Some compatibility checks and error message
  # Check whether the input A is a matrix
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("A should be an adjacency matrix with numeric entry of the graph.")
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
    stop("All the diagonal entries of A should be 0.")
  }

  #
  if(is.null(C) && is.null(K)){
    stop("Either block assignment or no of blocks should be provided.")
  }

  if (is.null(C)){
    # Getting block assignment for each node from adjacency matrix when block assignment is not provided
    C <- block_est(A, K)
  }

  if(is.null(K)){
    # Getting no of blocks from the block assignment vector
    K <- length(unique(C))
  }

  # Check all the elements of C and K is numeric
  if(!is.numeric(C) || !is.numeric(K)){
    stop("All the elements of C and K should be numeric.")
  }

  # Check whether all the elements of block assignment are from 1:K
  if (!all(C %in% 1:K)) {
    stop("C can only contain values from 1 to K.")
  }

  # Check whether there is at least one node from each of the block
  if(length(unique(C)) != K){
    stop("All the blocks should have atleast one node.")
  }

  # Getting dimension information
  n <- nrow(A)
  if (length(C) != n) {
    stop("The C should have same length as no of rows in A")
  }

  # Getting an igraph (an undirected and unweighted) object from the input adjacency matrix
  G_obs <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = NULL)

  # Calculate estimate of the parameter from observed graph
  # which will remain same after generating a new graph on same fiber
  p_mle <- get_mle(G_obs, C)

  # It will store GoF test statistic on graphs
  chi_seqR <- rep(0, numGraphs)
  G <- G_obs

  # Storing the first entry of chi_seq as test-stat on observed graph
  chi_seqR[1] <- round(graphchi(G, C, p_mle), 2)
  for (i in 2:numGraphs) {
    # Sampling a new graph
    G_current <- sample_a_move(C, G)

    # Computing GoF test statistic on new sampled graph
    chi_seqR[i] <- round(graphchi(G_current, C, p_mle), 2)
    G <- G_current
  }

  # pvalue i.e, proportion of sampled grpahs has larger GoF statistic than observed one
  pvalue <- mean(chi_seqR > chi_seqR[1])

  # Output:
  # chi_seq: sequence of chi square test statistics on the sampled graphs
  # pvalue: estimated p-value when true model is ERSBM
  return(list(statistic = chi_seqR, p.value = pvalue))
}
