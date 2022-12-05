#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Getting edge list in matrix
// adj - adjacency matrix of a graph
// E - no of edges in the graph
// n - no of nodes in the graph
// [[Rcpp::export]]
arma::mat get_edgelist(const arma::mat& adj, const int E, const int n){
  // Initialize edge_mat with all zeros
  arma::mat edge_mat(E, 2); edge_mat.zeros();

  // If there is an edge between two nodes, it will add those two nodes in edge_mat
  int e = 0;
  for (int row = 0; row < n; row ++) {
    for (int col = n - 1; col > row; col --) {
      if (adj(row, col) == 1){
        edge_mat(e, 0) = row + 1;
        edge_mat(e, 1) = col + 1;
        e = e + 1;
      }
    }
  }

  // Return output
  // edge_mat - a E x 2 matrix containing the edges of the graph
  return edge_mat;
}

// Getting edge list of a complement graph in matrix
// adj - adjacency matrix of a graph
// E - no of edges in the complement graph
// n - no of nodes in the graph
// [[Rcpp::export]]
arma::mat get_edgelist_comp(const arma::mat& adj, const int E, const int n){
  // Initialize edge_mat with all zeros
  arma::mat edge_mat(E, 2); edge_mat.zeros();

  // If there is no edge between two nodes i.e., there is an edge in the complement graph
  // it will add those two nodes in edge_mat
  int e = 0;
  for (int row = 0; row < n; row ++) {
    for (int col = n - 1; col > row; col --) {
      if (adj(row, col) == 0){
        edge_mat(e, 0) = row + 1;
        edge_mat(e, 1) = col + 1;
        e = e + 1;
      }
    }
  }

  // Return output
  // edge_mat - a E x 2 matrix containing the edges of the complement graph
  return edge_mat;
}
