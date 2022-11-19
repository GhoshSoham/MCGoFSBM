# The true model is ERSBM
# Chi-Square goodness of fit test statistic computation using the observed graph and MLE

#' Title
#'
#' @param G_obs igraph object which is an undirected graph and has no self loop
#' @param C numeric vector of size n of block assignment; from 1 to k
#' @param numGraphs number of graphs will be sampled; default value is 100
#'
#' @return chi_seq: sequence of chi square test statistics on the sampled graphs
#' @return pvalue: estimated p-value when true model is ERSBM
#' @export
#'
#' @examples
goftest <- function(G_obs, C, numGraphs = 100) {
  if (!inherits(G_obs, "igraph")) {
    stop("The input G_obs should be an igraph object")
  }

  if (igraph::is.weighted(G_obs)) {
    stop("The input graph G_obs shouldn't be a weighted one")
  }

  if (igraph::is_directed(G_obs)) {
    stop("The input graph G_obs shouldn't be a directed one")
  }

  if (igraph::any_loop(G_obs)) {
    stop("The graph G_obs shouldn't have any self loop")
  }

  n <- igraph::gorder(G_obs)
  if (length(C) != n) {
    stop("The C should have same length as no of vertices of graph")
  }

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
  return(list(chi_seq, pvalue))
}
