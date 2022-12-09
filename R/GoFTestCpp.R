#'
#' @title Monte Carlo Goodness of Fit tests for Stochastic Block Models
#'
#' @description `goftest_cpp` performs chi square goodness of fit test for network data under the model as ERSBM using  Rcpp
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
goftest_cpp <- function(A, K, C, numGraphs = 100) {
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

  # Check whether C or K is provided
  if(is.null(C) && is.null(K)){
    stop("Either block assignment or no of blocks should be provided.")
  }

  # If C is not provided, estimate C using value of K, the no of blocks
  if (is.null(C)){
    # Getting block assignment for each node from adjacency matrix when block assignment is not provided
    C <- block_est(A, K)
  }

  # If K is not provided, getting no of blocks from C
  if(is.null(K)){
    # Getting no of blocks from the block assignment vector
    K <- length(unique(C))
  }

  # Check all the elements of C and K are numeric
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

  # Check whether length of C is compatible with dimension of A
  if (length(C) != n) {
    stop("The C should have same length as no of rows in A")
  }

  # Calculating GoF statistic several times after sampling new graph on same fiber
  chi_seq_cpp <- graph_chain_on_fiber(A, C, n, K, numGraphs)
  chi_seq_cpp <- round(chi_seq_cpp[, 1], 2)

  # pvalue i.e, proportion of sampled grpahs has larger GoF statistic than observed one
  pvalue <- mean(chi_seq_cpp > chi_seq_cpp[1])

  # Output:
  # chi_seq: sequence of chi square test statistics on the sampled graphs
  # pvalue: estimated p-value when true model is ERSBM
  return(list(statistic = chi_seq_cpp, p.value = pvalue))
}
