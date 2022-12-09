#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Getting edge list in matrix
// adj - adjacency matrix of a graph
// E - no of edges in the graph
// n - no of nodes in the graph
// [[Rcpp::export]]
arma::umat get_edgelist_cpp(const arma::mat& adj, const arma::vec& C, const int E, const int n){
  // Initialize edge_mat with all zeros
  arma::umat edge_mat(E, 4); edge_mat.zeros();

  // If there is an edge between two nodes, it will add those two nodes in edge_mat
  int e = 0;
  for (int row = 0; row < n - 1; row ++) {
    for (int col = row + 1; col < n; col ++) {
      if (adj(row, col) == 1){
        edge_mat(e, 0) = row + 1;
        edge_mat(e, 1) = col + 1;
        edge_mat(e, 2) = C(row);
        edge_mat(e, 3) = C(col);
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
arma::umat get_edgelist_comp_cpp(const arma::mat& adj, const arma::vec& C, const int E, const int n){
  // Initialize edge_mat with all zeros
  arma::umat edge_mat(E, 4); edge_mat.zeros();

  // If there is no edge between two nodes i.e., there is an edge in the complement graph
  // it will add those two nodes in edge_mat
  int e = 0;
  for (int row = 0; row < n - 1; row ++) {
    for (int col = row + 1; col < n; col ++) {
      if (adj(row, col) == 0){
        edge_mat(e, 0) = row + 1;
        edge_mat(e, 1) = col + 1;
        edge_mat(e, 2) = C(row);
        edge_mat(e, 3) = C(col);
        e = e + 1;
      }
    }
  }

  // Return output
  // edge_mat - a E x 2 matrix containing the edges of the complement graph
  return edge_mat;
}

