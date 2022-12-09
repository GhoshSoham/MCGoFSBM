#'
#' @title Monte Carlo Goodness of Fit tests for Stochastic Block Models
#'
#' @description `goftest_cpp` performs chi square goodness of fit test for network data under the model as ERSBM using  Rcpp
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
goftest_cpp <- function(A, K, numGraphs = 100) {
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

  # Dimension information
  n <- nrow(A)

  # Getting block assignment for each node from adjacency matrix
  C_cpp <- block_est(A, K)

  # Calculating GoF statistic several times after sampling new graph on same fiber
  chi_seq_cpp <- graph_chain_on_fiber(A, C_cpp, n, K, numGraphs)
  chi_seq_cpp <- round(chi_seq_cpp[, 1], 2)

  # pvalue i.e, proportion of sampled grpahs has larger GoF statistic than observed one
  pvalue <- mean(chi_seq_cpp > chi_seq_cpp[1])

  # Output:
  # chi_seq: sequence of chi square test statistics on the sampled graphs
  # pvalue: estimated p-value when true model is ERSBM
  return(list(statistic = chi_seq_cpp, p.value = pvalue))
}
