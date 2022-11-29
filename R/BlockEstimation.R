# A a n by n binary symmetric adjacency matrix representing a undirected graph where n is the no nodes in the graph
# K a numeric scalar representing no of blocks
#
# cluster: a vector of size n representing block assignment for each node; values are 1 to K i.e, no of cluster

block_est <- function (A, K)
{
  SVD <- svd(A, nu = K, nv = K)

  km <- stats::kmeans(SVD$v[, 1:K], centers = K, nstart = 30, iter.max = 30)
  return(cluster = km$cluster)
}
