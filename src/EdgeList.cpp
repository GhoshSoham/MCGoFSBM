#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Getting edge list in matrix
// adj - adjacency matrix of a graph
// E - no of edges in the graph
// n - no of nodes in the graph
// [[Rcpp::export]]
arma::mat get_edgelist(const arma::mat& adj, const int E, const int n){
  // Initialize beta_mat and beta_start with all zeros
  arma::mat edge_mat(E, 2); edge_mat.zeros();

  int e = 0;
  for (int row = 0; row<n; row++) {
    for (int col = n-1; col>row; col--) {
      if (adj(row, col) == 1){
        edge_mat(e, 0) = row + 1;
        edge_mat(e, 1) = col + 1;
        e = e+1;
      }
    }
  }

  // Return output
  // edge_mat - a E x 2 matrix containing the edges of the graph
  return edge_mat;
}
