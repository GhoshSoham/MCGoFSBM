# The true model is ERSBM
# Testing whether the observed graph follows ERSBM model using chi square goodness of fit test
# for network data

#'
#' @title Monte Carlo Goodness of Fit tests for Stochastic Block Models
#'
#' @description `goftest` performs chi square goodness of fit test for network data considering the model as ERSBM
#'
#' @param A a n by n binary symmetric adjacency matrix representing a undirected graph where n is the no nodes in the graph
#' @param K a positive integer scalar representing no of blocks; K>1
#' @param C an positive integer vector of size n of block assignment of each node; from 1 to K (no of blocks)
#' @param numGraphs number of graphs will be sampled; default value is 100
#'
#' @return A list with the elements
#' \item{statistic}{the values of the chi-square test statistics on each sampled graph}
#' \item{p.value}{the p-value for the test.}
#'
#' @export
#'
#' @examples
#' RNGkind(sample.kind = "Rounding")
#' set.seed(1729)
#'
#' # We model a network with 3 even classes.
#' n1 <- 50
#' n2 <- 50
#' n3 <- 50
#'
#' # Generating block assignment for each of the nodes
#' n <- n1 + n2 + n3
#' class <- rep(c(1, 2, 3), c(n1, n2, n3))
#'
#' # Generating the adjacency matrix of the network
#' # Generate the matrix of connection probability.
#' cmat <- matrix(
#'   c(
#'     30, 0.05, 0.05,
#'     0.05, 30, 0.05,
#'     0.05, 0.05, 30
#'   ),
#'   ncol = 3,
#'   byrow = TRUE
#' )
#' pmat <- cmat / n
#'
#' # Creating the n x n adjacency matrix.
#' adj <- matrix(0, n, n)
#' for (i in 2:n) {
#'   for (j in 1:(i - 1)) {
#'     p <- pmat[class[i], class[j]] # We find the probability of connection with the weights.
#'     adj[i, j] <- rbinom(1, 1, p) # We include the edge with probability p.
#'   }
#' }
#'
#' adjsymm <- adj + t(adj)
#'
#' # When class assignment is known
#' out <- goftest(adjsymm, C = class, numGraphs = 100)
#'
#' chi_sq_seq <- out$statistic
#' pvalue <- out$p.value
#' print(pvalue)
#'
#' # Plotting histogram of the sequence of the test statistics
#' hist(chi_sq_seq, 20, xlab = "chi-square test statistics", main = NULL)
#' abline(v = chi_sq_seq[1], col = "red", lwd = 5) # adding test statistic on the observed network
#'
#' @references
#' Karwa et al. (2016) Monte Carlo goodness-of-fit tests for degree corrected and related stochastic blockmodels,
#' \doi{10.48550/ARXIV.1612.06040}
#'
#' @references
#' Qin, T., & Rohe, K. (2013). Regularized spectral clustering under the degree-corrected stochastic blockmodel,
#' \emph{Advances in neural information processing systems},
#' \strong{26}.
#'
#' @references
#' Lei, J., & Rinaldo, A. (2015). Consistency of spectral clustering in stochastic block models,
#' \emph{The Annals of Statistics},
#' \strong{43(1)}, 215-237.
#' \doi{10.1214/14-aos1274}
#'
#' @references
#' Li T, Levina E, Zhu J (2021). randnet: Random Network Model Estimation, Selection and Parameter Tuning.
#' \emph{R packageversion 0.3},
#' <https://CRAN.R-project.org/package=randnet>.
#'

goftest <- function(A, K = NULL, C = NULL, numGraphs = 100) {
  # Some compatibility checks and error message
  # Check whether the input A is a matrix
  if (!is.matrix(A) | !is.numeric(A)) {
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

  # Check whether numGraphs is numeric and positive
  if (!is.numeric(numGraphs) || numGraphs <= 0 || numGraphs %% 1 != 0) {
    stop("numGraphs, number of graphs to sample should be a positive natural number.")
  }

  # Check whether C or K is provided
  if (is.null(C) & is.null(K)) {
    stop("Either block assignment or no of blocks should be provided.")
  }

  # If K is not provided, getting no of blocks from C
  if (is.null(K)) {
    # Getting no of blocks from the block assignment vector
    K <- length(unique(C))
  }

  # Check whether K is numeric and positive integer
  if (!is.numeric(K) | K %% 1 != 0) {
    stop("No of blocks K should be a positive natural number.")
  }

  # No of blocks K should be at least 2
  if (K < 2) {
    stop("No of blocks K should be at least 2.")
  }

  # If C is not provided, estimate C using value of K, the no of blocks
  if (is.null(C)) {
    # Getting block assignment for each node from adjacency matrix when block assignment is not provided
    C <- block_est(A, K)
  }

  # Check whether C is numeric
  if (!is.numeric(C)) {
    stop("All the elements of C should be numeric.")
  }

  # Check whether all the elements of block assignment are from 1:K
  if (!all(C %in% 1:K)) {
    stop("C can only contain values from 1 to K (no of blocks).")
  }

  # Check whether there is at least one node from each of the block
  if (length(unique(C)) != K) {
    stop("All the blocks should have atleast one node.")
  }

  # Getting dimension information
  n <- nrow(A)

  # Check whether length of C is compatible with dimension of A
  if (length(C) != n) {
    stop("The C should have same length as no of rows in A.")
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
